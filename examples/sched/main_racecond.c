#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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

// This is the same example as in main_task_model.c, but with 
// an artificial race condition. Apart from the control threads
// we start another one who fills all columns with zeros while 
// the controllers are at work.

////////////////////////////////////////////////////////////////////////////////////////
// some constants for the example                                                     //
////////////////////////////////////////////////////////////////////////////////////////

// random numbers are picked between 0 and RND_MAX
#define RND_MAX 10000
// number of columns to be filled (=number of control threads)
#define NUM_TASKS 4
// we perform three different operations, rndX, incX and divX
#define NUM_JOBTYPES 3 
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
  
  // these are all blocking function calls right now
  v[0]=0;// next one will be randomized
  for (i=1; i<n;i++)
    {
    v[i]=v[i-1];
    if (v[i-1]==0 || v[i-1]==1)
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

// create an artificial race condition by writing in the memory
// that is being filled by the fill_vector function, thus messing
// up the results. When the executing threads encounter the false
// zeros, they will randomize and restart the sequence falsely. If
// the zero is placed after the entry was generated, the sequence 
// is wrong. If the zero is written before the entry is generated,
// nothing happens.
void* init_vector(void* v)
  {
  int i,j;
  int* vectors = (int*)v;
  for (i=0; i<NDIM;i++)
    {
    for (j=0;j<NUM_TASKS;j++)
      {
      // delay for up to 100 microns
      usleep((useconds_t)(rand()/RAND_MAX*100));
      vectors[j*NDIM+i]=0;
      }
    }
  fprintf(OUT,"vector initialized\n");
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
#ifndef PIN_CONTROL_THREADS
  pthread_t controlThread[NUM_TASKS];
#else
  ghost_task_t *controlTask[NUM_TASKS];
#endif
  pthread_t initThread;
  
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

#ifndef PIN_CONTROL_THREADS

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

// this thread sets random entries in the vector to 0, which it should do
// (if at all) before the other starts. The lack of synchronization here 
// creates a race condition.
pthread_create(&initThread, NULL, &init_vector, (void *)vectors);


for (i=0;i<NUM_TASKS;i++)
  {
  mainArg[i].n=NDIM;
  mainArg[i].col=i;
  mainArg[i].v=vectors+i*NDIM;
  mainArg[i].taskBuf=taskBuf;
  mainArg[i].RNDX=op_RNDX;
  mainArg[i].INCX=op_INCX;
  mainArg[i].DIVX=op_DIVX;

#ifdef PIN_CONTROL_THREADS
  controlTask[i] = ghost_task_init(1, 0, &fill_vector, 
        (void*)(&mainArg[i]),GHOST_TASK_DEFAULT);
  ghost_task_add(controlTask[i]);
#else
  pthread_create(&controlThread[i], NULL, &fill_vector, (void *)&mainArg[i]);
#endif
  }

#ifdef PIN_CONTROL_THREADS
ghost_task_waitall();
#else

pthread_join(initThread,NULL);

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
