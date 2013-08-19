#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <pthread.h>
#include <time.h>

#include "phist_taskbuf.h"

#include "phist_kernels.h"
#include "phist_macros.h"

#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#include "ghost_taskq.h"
#else
#error "this driver makes no sense without ghost"
#endif

//                                                              
// this example driver performs the following algorithm.        
// given a number of columns of integer vectors, fill each      
// vector with the sequence                                     
//                                                              
// V[n+1]=random() if n=-1 or V[n]=0 or V[n]=1                  
// V[n+1] = V[n]+1 if V[n] is odd                               
// V[n+1] = V[n]/2 if n is even                                 
//                                                              
// Each vector is treated by a different 'control thread', the  
// actual operations (rndX(), incX and divX) are put into a     
// task buffer and executed in bulks whenever each control      
// thread has announced the type of operation it needs.         
//                                                              

////////////////////////////////////////////////////////////////////////////////////////
// some constants for the example                                                     //
////////////////////////////////////////////////////////////////////////////////////////

// random numbers are picked between 0 and RND_MAX
#define RND_MAX 10000
// number of columns to be filled (=number of control threads)
#define NUM_TASKS 4
// local length of each column
#define NDIM 40

#define OUT stderr

////////////////////////////////////////////////////////////////////////////////////////
// algorithm-specific functionality                                                   //
////////////////////////////////////////////////////////////////////////////////////////

// fill  all entries in an array with random integers
static void *rndX(argList_t* args)
  {
  int i;

//  fprintf(OUT,"NUM THREADS: %d\n",arg->nThreads);
  omp_set_num_threads(args->nthreads);

#pragma omp parallel
  {
#pragma omp single
  fprintf(OUT,"rndX executed by %d threads on %d values\n",omp_get_num_threads(),args->n);  
#pragma omp for
  for (i=0;i<args->n;i++)
    {
    int* val = (int*)args->arg[i];
    *val=(int)((rand()/(double)RAND_MAX)*RND_MAX);
    }
  }
  return NULL;
  }

// increment all entries in an array by one
static void *incX(argList_t* args)
  {
  int i;

//  fprintf(OUT,"NUM THREADS: %d\n",arg->nThreads);
  omp_set_num_threads(args->nthreads);

#pragma omp parallel
  {
#pragma omp single
  fprintf(OUT,"incX executed by %d threads on %d values\n",
        omp_get_num_threads(),args->n);  
#pragma omp for
  for (i=0;i<args->n;i++)
    {
    int* val = (int*)args->arg[i];
    *val+=1;
    }
  }
  return NULL;
  }


// divide all entries in an array by two
static void *divX(argList_t* args)
  {
  int i;

//  fprintf(OUT,"NUM THREADS: %d\n",arg->nThreads);
  omp_set_num_threads(args->nthreads);

#pragma omp parallel
  {
#pragma omp single
  fprintf(OUT,"divX executed by %d threads on %d values\n", 
        omp_get_num_threads(),args->n);  
#pragma omp for
  for (i=0;i<args->n;i++)
    {
    int* val = (int*)args->arg[i];
    *val/=2;
    }
  }
  return NULL;
  }

////////////////////////////////////////////////////////////////////////////////////////////
// implementation of the algorithm                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

// to conform with the ghost queue, we pass all function arguments to
// our "algorithm" fill_vector in a single void pointer:
typedef struct {
int n;
int col;
int *v;
int RNDX; // key for calling rndX via the buffer
int INCX; // key for calling incX via the buffer
int DIVX; // key for calling divX via the buffer
taskBuf_t *taskBuf; // task buffer to be used, must be initialized
                    // to handle the above three operations
int ierr; // return code
} mainArg_t;

// this is a long-running job (an "algorithm") executed by a control thread.
// Inside, work that has to be done is put in a task buffer and than bundled
// between control threads and put into the ghost queue.
void* fill_vector(void* v_mainArg)
  {
  int i,n,col;
  int* v;
  int* ierr;
  taskBuf_t* taskBuf;
  mainArg_t *mainArg = (mainArg_t*)v_mainArg;

  n = mainArg->n;
  col = mainArg->col;
  v = mainArg->v;
  taskBuf=mainArg->taskBuf;
  ierr = &mainArg->ierr;
  
  fprintf(OUT,"thread %ul: taskBuf at %p, finished_cv at %p\n",pthread_self(),
        taskBuf, &taskBuf->finished_cv);
  
  // these are all blocking function calls right now
  v[0]=0;// next one will be randomized
  for (i=1; i<n;i++)
    {
    v[i]=v[i-1];
    if (v[i-1]==42)
      {
      PHIST_CHK_IERR(taskBuf_renounce(taskBuf,col,ierr),*ierr);
      return NULL;
      }
    else if (v[i-1]==0 || v[i-1]==1)
      {
      PHIST_CHK_IERR(taskBuf_add(taskBuf,&v[i],col,mainArg->RNDX,ierr),*ierr);
      }
    else if (v[i-1]%2==0)
      {
      PHIST_CHK_IERR(taskBuf_add(taskBuf,&v[i],col,mainArg->DIVX,ierr),*ierr);
      }
    else
      {
      PHIST_CHK_IERR(taskBuf_add(taskBuf,&v[i],col,mainArg->INCX,ierr),*ierr);
      }
#ifndef SYNC_WAIT
    PHIST_CHK_IERR(taskBuf_wait(taskBuf,col,ierr),*ierr);
#endif
    }
  return NULL;
  }

  
int main(int argc, char** argv)
  {
  int ierr;
  int i,j;
  int nworkers;
  int op_RNDX, op_INCX, op_DIVX;
  taskBuf_t *taskBuf;
  mainArg_t mainArg[NUM_TASKS];
#ifdef DONT_PIN_CONTROL_THREADS
  pthread_t controlThread[NUM_TASKS];
#else
  ghost_task_t *controlTask[NUM_TASKS];
#endif
  
  int* vectors = (int*)malloc(NUM_TASKS*NDIM*sizeof(int));

  // we don't really use the kernel lib in this example, but
  // out of habit we initialize it:
  PHIST_CHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);
// avoid double init
#ifndef PHIST_KERNEL_LIB_GHOST  
  ghost_init(argc,argv);
#endif

  // initialize C random number generator
  srand(time(NULL));

#ifdef DONT_PIN_CONTROL_THREADS

nworkers=ghost_thpool->nThreads;

#else
// we would actually like to run the control threads without pinning
// and allow having all cores for doing the work, but if we do that 
// in the current ghost implementation, the jobs in the ghost queue 
// never get executed because NUM_TASKS cores are reserved for the     
// control threads.
nworkers=ghost_thpool->nThreads - NUM_TASKS;

#endif

PHIST_CHK_IERR(taskBuf_create(&taskBuf,NUM_TASKS,&ierr),ierr);

// add the basic operations for our algorithm to the buffer
PHIST_CHK_IERR(taskBuf_add_op(taskBuf, &rndX, nworkers,&op_RNDX, &ierr),ierr);
PHIST_CHK_IERR(taskBuf_add_op(taskBuf, &incX, nworkers,&op_INCX, &ierr),ierr);
PHIST_CHK_IERR(taskBuf_add_op(taskBuf, &divX, nworkers,&op_DIVX, &ierr),ierr);

for (i=0;i<NUM_TASKS;i++)
  {
  mainArg[i].n=NDIM;
  mainArg[i].col=i;
  mainArg[i].v=vectors+i*NDIM;
  mainArg[i].taskBuf=taskBuf;
  mainArg[i].RNDX=op_RNDX;
  mainArg[i].INCX=op_INCX;
  mainArg[i].DIVX=op_DIVX;
#ifndef DONT_PIN_CONTROL_THREADS
  controlTask[i] = ghost_task_init(1, 0, &fill_vector, 
        (void*)(&mainArg[i]),GHOST_TASK_DEFAULT);
  ghost_task_add(controlTask[i]);
#else
  pthread_create(&controlThread[i], NULL, &fill_vector, (void *)&mainArg[i]);
#endif
  }

#ifndef DONT_PIN_CONTROL_THREADS
ghost_task_waitall();
#else
for (i=0;i<NUM_TASKS;i++)
  {
  pthread_join(controlThread[i],NULL);
  }
#endif

for (i=0;i<NUM_TASKS;i++)
  {
  if (mainArg[i].ierr)
    {
    fprintf(stderr,"error flag %d returned by thread %d\n",mainArg[i].ierr,i);
    ierr=mainArg[i].ierr; // the last non-zero return code will be returned
    }
  }
// single thread prints the results
for (i=0;i<NDIM;i++)
  {
  for (j=0;j<NUM_TASKS;j++)
    {
    fprintf(OUT,"\t%d",vectors[j*NDIM+i]);
    }
  fprintf(OUT,"\n");
  }
  
  PHIST_CHK_IERR(phist_kernels_finalize(&ierr),ierr);
  // avoid duplicate finish call
#ifndef PHIST_KERNEL_LIB_GHOST
  // finalize ghost queue
  ghost_finish();
#endif  
  }
