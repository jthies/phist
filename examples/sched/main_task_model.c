#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <pthread.h>
#include <time.h>

#include "phist_taskbuf.h"

#include "phist_macros.h"

#include "ghost.h"
#include "ghost_taskq.h"

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
// operations incX and divX are put into a                      
// task buffer and executed in bulks whenever each control      
// thread has announced the type of operation it needs.         
// The operation rndX is not buffered but enqueued on the number
// of cores each control thread gets for local operations.      
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

////////////////////////////////////////////////////////////////////////////////////////
// algorithm-specific functionality                                                   //
////////////////////////////////////////////////////////////////////////////////////////

// fill  all entries in an array with random integers
static void rndX(int* arg)
  {

#pragma omp parallel
  {
#pragma omp single
    {
    PHIST_OUT(0,"rndX executed by %d threads on one value\n",omp_get_num_threads());  
    *arg=(int)((rand()/(double)RAND_MAX)*RND_MAX);
    }
  }
  return;
  }

// increment all entries in an array by one
static void incX(argList_t* args)
  {
  int i;

//  PHIST_OUT(0,"NUM THREADS: %d\n",arg->nThreads);

#pragma omp parallel
  {
#pragma omp single
  PHIST_OUT(0,"incX executed by %d threads on %d values\n",
        omp_get_num_threads(),args->n);  
#pragma omp for
  for (i=0;i<args->n;i++)
    {
    *((int*)args->out_arg[i])=*((const int*)args->in_arg[i])+1;
    }
  }
  return;
  }


// divide all entries in an array by two
static void divX(argList_t* args)
  {
  int i;

//  PHIST_OUT(0,"NUM THREADS: %d\n",arg->nThreads);

#pragma omp parallel
  {
#pragma omp single
  PHIST_OUT(0,"divX executed by %d threads on %d values\n", 
        omp_get_num_threads(),args->n);  
#pragma omp for
  for (i=0;i<args->n;i++)
    {
    *((int*)args->out_arg[i])=*((const int*)args->in_arg[i])/2;
    }
  }
  return;
  }

////////////////////////////////////////////////////////////////////////////////////////////
// implementation of the algorithm                                                        //
///////////////////
/////////////////////////////////////////////////////////////////////////

// to conform with the ghost queue, we pass all function arguments to
// our "algorithm" fill_vector in a single void pointer:
typedef struct {

int nworkers;
int n;
int col;
int *v;
int RNDX; // key for calling rndX via the buffer
int INCX; // key for calling incX via the buffer
int DIVX; // key for calling divX via the buffer
taskBuf_t *taskBuf; // task buffer to be used, must be initialized
                    // to handle the above three operations
int ierr; // return code
} colArg_t;

// to conform with the ghost queue, we pass all function arguments to
// our "algorithm" fill_all in a single void pointer:
typedef struct {

int nworkers;
int n; // length of the vectors
int nv; // number of vectors/columns
int *v; // pointer to the first vector

int ierr; // return code
} mainArg_t;

// forward declaration because fill_all needs fill_vector
void fill_vector(colArg_t* v_colArg);


void* fill_all(mainArg_t* arg)
  {
  colArg_t colArg[arg->nv];
  ghost_task_t *colTask[arg->nv];
  int nteam; // number of threads per column that can do local operations
  int op_INCX; // key for calling incX via the buffer
  int op_DIVX; // key for calling divX via the buffer
  taskBuf_t *taskBuf; // task buffer to be used, must be initialized
                      // to handle the above three operations
  int *ierr = &arg->ierr;
  int i;
                 
  nteam = (int) ( (double)arg->nworkers / (double)arg->nv );                    
  
  if (nteam<1)
    {
    PHIST_OUT(0,"number of threads (%d) smaller than number of tasks (%d)\n",
        arg->nworkers, arg->nv);
    arg->ierr=-1; 
    return;
    }
  // first, setup the task buffer
  PHIST_CHK_IERR(taskBuf_create(&taskBuf,arg->nv,ierr),*ierr);

  // add the basic operations for our algorithm to the buffer
  PHIST_CHK_IERR(taskBuf_add_op(taskBuf, (void*(*)(argList_t*))&incX, NULL, 
        arg->nworkers, &op_INCX, ierr),*ierr);
  PHIST_CHK_IERR(taskBuf_add_op(taskBuf, (void*(*)(argList_t*))&divX, NULL, 
        arg->nworkers, &op_DIVX, ierr),*ierr);

  // enqueue the tasks of filling the columns
  for (i=0;i<arg->nv;i++)
    {
    colArg[i].nworkers=nteam;
    colArg[i].n=arg->n;
    colArg[i].col=i;
    colArg[i].v=arg->v+i*arg->n;
    colArg[i].taskBuf=taskBuf;
    colArg[i].INCX=op_INCX;
    colArg[i].DIVX=op_DIVX;

    colTask[i] = ghost_task_init(nteam, GHOST_TASK_LD_ANY, (void*(*)(void*))&fill_vector, 
        (void*)(&colArg[i]),GHOST_TASK_DEFAULT|GHOST_TASK_NO_PIN);
    ghost_task_add(colTask[i]);
    }
  // now we have started all the teams of threads that will work
  // on a column and request global operations from the buffer. 
  // The teams themselves kick off the buffer ops whenever it is
  // full. As this master thread has currently control over all 
  // the available cores, it now has to wait explicitly for the 
  // started control threads so that it's resources can be used 
  // by them (otherwise they can't start running).
  for (i=0;i<arg->nv;i++)
    {
    PHIST_OUT(0,"master thread %lu waiting for col %d\n",pthread_self(),i);
    ghost_task_wait(colTask[i]);
    ghost_task_destroy(colTask[i]);
    }

  PHIST_OUT(0,"master thread %lu has finished, so all controllers have renounced the buffer",
        pthread_self());
  PHIST_CHK_IERR(taskBuf_delete(taskBuf,ierr),*ierr);
  return;
  }

// this is a long-running job (an "algorithm") executed by a control thread.
// Inside, work that has to be done is put in a task buffer and than bundled
// between control threads and put into the ghost queue.
void fill_vector(colArg_t* colArg)
  {
  int i,n,col;
  int* v;
  int* ierr;
  taskBuf_t* taskBuf;
  ghost_task_t *rndX_task;
  int no_wait;

  n = colArg->n;
  col = colArg->col;
  v = colArg->v;
  taskBuf=colArg->taskBuf;
  ierr = &colArg->ierr;
  
#pragma omp parallel
#pragma omp single
PHIST_OUT(0,"fill_vector on col %d, team size %d, task buffer at %p\n",
        col, omp_get_num_threads(), taskBuf);

//  PHIST_OUT(0,"thread %ul: taskBuf at %p, finished_cv at %p\n",pthread_self(),
//        taskBuf, &taskBuf->finished_cv);

  // these are all blocking function calls right now
  v[0]=0;// next one will be randomized
  
  for (i=1; i<n;i++)
    {
    rndX_task=NULL;
    no_wait=0;
    if (v[i-1]==42)
      {
      PHIST_CHK_IERR(taskBuf_renounce(taskBuf,col,ierr),*ierr);
      return;
      }
    else if (v[i-1]==0 || v[i-1]==1)
      {
      // this is not done via the buffer but using the own team of OpenMP threads

      // create a ghost task for the rndX operation
      // TODO - keep track of locality domains (LD) and do this correctly
      // TODO - is it safe to change the arg but keep the task?
      // NOTE: with the current implementation we could even run the rndX
      //       function locally and skip the queuing completely, but that
      //       may keep others waiting. In the unlikely situation that one
      //       thread just performs rndX all the time it will finish the
      //       whole column before anyone else can start doing something.
      //       Using the queue here allows the buffered tasks to get in 
      //       between the local operations.
      
      
      // the variant with the queue doesn't work yet, so we try executing rndX
      // locally first
      rndX(&v[i]);
      no_wait=1;
      /*
      rndX_task = ghost_task_init(colArg->nworkers, GHOST_TASK_LD_ANY,
        (void*(*)(void*))rndX, (void*)(&v[i]), GHOST_TASK_DEFAULT);
      ghost_task_add(rndX_task);
      */
 
      // we have two options here: tell the buffer to do this round without us,
      // or keep all other controllers waiting while we do the rndX. The choice
      // depends on the priority of the rndX operation.
      
      // TODO - right now there is no mechanism to prevent the same thread from putting
      // in another task on the column, in which case the countdown is decremented twice
      // by the same controller
//      PHIST_CHK_IERR(taskBuf_add(taskBuf,NULL,NULL,col,PHIST_OP_NULL,ierr),*ierr);
      }
    else if (v[i-1]%2==0)
      {
      PHIST_CHK_IERR(taskBuf_add(taskBuf,&v[i-1],&v[i],col,colArg->DIVX,ierr),*ierr);
      }
    else
      {
      PHIST_CHK_IERR(taskBuf_add(taskBuf,&v[i-1],&v[i],col,colArg->INCX,ierr),*ierr);
      }
    // if we did not put in a task this iteration, this function
    // will return immediately.
    if (no_wait!=1)
      {
      PHIST_CHK_IERR(taskBuf_wait(taskBuf,col,ierr),*ierr);
      }

    if (rndX_task!=NULL)
      {
      // wait for the ghost tasks started by ourself (rndX)
      ghost_task_wait(rndX_task);
      }
    }
  return;
  }

  
int main(int argc, char** argv)
  {
  int ierr;
  int i,j;
  int nworkers;
  
  mainArg_t mainArg;
  ghost_task_t *mainTask;
  
  int* vectors = (int*)malloc(NUM_TASKS*NDIM*sizeof(int));

  PHIST_OUT(0,"initialize ghost");
  ghost_init(argc,argv);
  ghost_thpool_init(GHOST_THPOOL_NTHREADS_FULLNODE,GHOST_THPOOL_FTHREAD_DEFAULT,GHOST_THPOOL_LEVELS_FULLSMT);
  ghost_taskq_init();

  // initialize C random number generator
  srand(time(NULL));

nworkers=ghost_thpool->nThreads;

  mainTask = ghost_task_init(nworkers, GHOST_TASK_LD_ANY, (void*(*)(void*))&fill_all, 
        (void*)&mainArg,GHOST_TASK_DEFAULT|GHOST_TASK_NO_PIN);
  mainArg.nworkers=nworkers;
  mainArg.nv = NUM_TASKS;
  mainArg.n  = NDIM;
  mainArg.v = vectors;
  ghost_task_add(mainTask);

ghost_task_wait(mainTask);

  ghost_task_destroy(mainTask);  

if (mainArg.ierr)
  {
  fprintf(stderr,"error flag %d returned by master thread\n",mainArg.ierr);
  ierr=mainArg.ierr; // the last non-zero return code will be returned
  }
// single thread prints the results
for (i=0;i<NDIM;i++)
  {
  for (j=0;j<NUM_TASKS;j++)
    {
    fprintf(stdout,"\t%d",vectors[j*NDIM+i]);
    }
  fprintf(stdout,"\n");
  }

  free(vectors);


  // finalize ghost queue
  ghost_finish();
  }
