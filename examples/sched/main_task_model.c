#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <pthread.h>

#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#include "ghost_taskq.h"
#else
#error "this driver makes no sense without ghost"
#endif

#define DONT_PIN_CONTROL_THREADS

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

// no operation
#define JOBTYPE_NOOP -1
// x=random()
#define JOBTYPE_RNDX 0
// x=x+1
#define JOBTYPE_INCX 1
// x=x/2
#define JOBTYPE_DIVX 2

static const int ntasks=4; // number of columns to be filled (=number of control threads)
static const int njobTypes=3; // we perform three different operations, rndX, incX and divX
static const int ndim=25; // local length of each column

static int nworkers;
static FILE* out;

////////////////////////////////////////////////////////////////////////////////////////
// algorithm-specific functionality                                                   //
////////////////////////////////////////////////////////////////////////////////////////

// argument for rndX, incX and divX functions operating on an array of int pointers.
typedef struct
  {
  int nThreads;
  int n;
  int **val;
  } argList_t;

void argList_create(argList_t** args)
  {
  *args = (argList_t*)malloc(sizeof(argList_t));
  (*args)->nThreads=nworkers;
  (*args)->n=0;
  (*args)->val=(int**)malloc(ntasks*sizeof(int*));
  }
  
void argList_delete(argList_t* args)
  {
  free(args->val);
  free(args);
  }

// fill  all entries in an array with random integers
static void *rndX(void* v_arg)
  {
  int i;
  argList_t* arg = (argList_t*)v_arg;

  fprintf(out,"NUM WORKERS: %d\n",nworkers);
  fprintf(out,"NUM THREADS: %d\n",arg->nThreads);
  omp_set_num_threads(arg->nThreads);

#pragma omp parallel
  {
#pragma omp single
  fprintf(out,"rndX executed by %d threads on %d values\n",omp_get_num_threads(),arg->n);  
#pragma omp for
  for (i=0;i<arg->n;i++)
    *(arg->val[i])=(int)((rand()/(double)RAND_MAX)*1000);
  }
  return NULL;
  }

// increment all entries in an array by one
static void *incX(void* v_arg)
  {
  int i;
  argList_t* arg = (argList_t*)v_arg;

  fprintf(out,"NUM WORKERS: %d\n",nworkers);
  fprintf(out,"NUM THREADS: %d\n",arg->nThreads);
  omp_set_num_threads(arg->nThreads);

#pragma omp parallel
  {
#pragma omp single
  fprintf(out,"incX executed by %d threads on %d values\n",
        omp_get_num_threads(),arg->n);  
#pragma omp for
  for (i=0;i<arg->n;i++)
    *(arg->val[i])+=1;
  }
  return NULL;
  }


// divide all entries in an array by two
static void *divX(void* v_arg)
  {
  int i;
  argList_t* arg = (argList_t*)v_arg;

  fprintf(out,"NUM WORKERS: %d\n",nworkers);
  fprintf(out,"NUM THREADS: %d\n",arg->nThreads);
  omp_set_num_threads(arg->nThreads);

#pragma omp parallel
  {
#pragma omp single
  fprintf(out,"divX executed by %d threads on %d values\n", 
        omp_get_num_threads(),arg->n);  
#pragma omp for
  for (i=0;i<arg->n;i++)
    *(arg->val[i])/=2;
  }
  return NULL;
  }

////////////////////////////////////////////////////////////////////////////////////////
// task buffer prototype                                                              //
////////////////////////////////////////////////////////////////////////////////////////

// a buffer object, when each control task has set its jobType flag,
// the rndX, incX and divX jobs are bundled and enqueued in the ghost queue.
typedef struct {

int njobs; // total number of jobs (operating on a single int value)
int countdown; // when this reaches 0, the jobs are launched
pthread_mutex_t lock_mx; // for controlling access to the shared object and the condition 
                         // variable finished_cv.
pthread_cond_t finished_cv; // for waiting for the buffer to be emptied

int *jobType;
int **jobArg;

ghost_task_t **ghost_task;
argList_t **ghost_args;

} taskBuf_t;

void taskBuf_create(taskBuf_t** buf, int num_tasks, int num_workers, int* ierr);
void taskBuf_delete(taskBuf_t* buf, int* ierr);
void taskBuf_exec(taskBuf_t* buf, int* arg, int task_id, int taskFlag);
void taskBuf_group_and_run(taskBuf_t* buf);

// create a new task buffer object
void taskBuf_create(taskBuf_t** buf, int num_tasks, 
        int num_workers, int* ierr)
  {
  int i;
  *ierr=0;
  *buf = (taskBuf_t*)malloc(sizeof(taskBuf_t));
  
  (*buf)->njobs=num_tasks;
  (*buf)->countdown=num_tasks;
  (*buf)->jobType=(int*)malloc(num_tasks*sizeof(int));
  (*buf)->jobArg=(int**)malloc(num_tasks*sizeof(int*));

  pthread_mutex_init(&((*buf)->lock_mx),NULL);
  pthread_cond_init(&(*buf)->finished_cv,NULL);

  (*buf)->ghost_args=(argList_t**)malloc(njobTypes*sizeof(argList_t*));
  (*buf)->ghost_task=(ghost_task_t**)malloc(njobTypes*sizeof(ghost_task_t*));

  for (i=0;i<njobTypes;i++)
    {
    argList_create(&(*buf)->ghost_args[i]);
    (*buf)->ghost_args[i]->nThreads = num_workers;
    }
#ifdef DONT_PIN_CONTROL_THREADS
  (*buf)->ghost_task[JOBTYPE_RNDX] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &rndX, 
        (void*)((*buf)->ghost_args[JOBTYPE_RNDX]),GHOST_TASK_DEFAULT);

  (*buf)->ghost_task[JOBTYPE_INCX] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &incX, 
        (void*)((*buf)->ghost_args[JOBTYPE_INCX]),GHOST_TASK_DEFAULT);

  (*buf)->ghost_task[JOBTYPE_DIVX] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &divX, 
        (void*)((*buf)->ghost_args[JOBTYPE_DIVX]),GHOST_TASK_DEFAULT);
#else
  (*buf)->ghost_task[JOBTYPE_RNDX] = ghost_task_init(num_workers, 0, &rndX, 
        (void*)((*buf)->ghost_args[JOBTYPE_RNDX]),GHOST_TASK_DEFAULT);

  (*buf)->ghost_task[JOBTYPE_INCX] = ghost_task_init(num_workers, 0, &incX, 
        (void*)((*buf)->ghost_args[JOBTYPE_INCX]),GHOST_TASK_DEFAULT);

  (*buf)->ghost_task[JOBTYPE_DIVX] = ghost_task_init(num_workers, 0, &divX, 
        (void*)((*buf)->ghost_args[JOBTYPE_DIVX]),GHOST_TASK_DEFAULT);
#endif
  return;
  }

// create a new task buffer object
void taskBuf_delete(taskBuf_t* buf, int* ierr)
  {
  int i;
  *ierr=0;
  //TODO - maybe we should force execution/wait here
  free(buf->jobType);
  free(buf->jobArg);
  for (i=0;i<njobTypes;i++)
    {
    ghost_task_destroy(buf->ghost_task[i]);
    argList_delete(buf->ghost_args[i]);
    }
  pthread_mutex_destroy(&buf->lock_mx);
  pthread_cond_destroy(&buf->finished_cv);
  free(buf->ghost_task);
  free(buf->ghost_args);
  free(buf);
  return;
  }

// a job is put into the task buffer and executed once the
// buffer is full. This function returns only when the task 
// is finished.
void taskBuf_exec(taskBuf_t* buf, int* arg, int task_id, int taskFlag)
  {
  // lock the buffer while putting in jobs and executing them
  pthread_mutex_lock(&buf->lock_mx);
  buf->jobType[task_id]=taskFlag;
  buf->jobArg[task_id]=arg;
  buf->countdown--;
  fprintf(out,"control thread %lu exec job %d on column %d, countdown=%d\n",
        pthread_self(), taskFlag, task_id, buf->countdown);
  if (buf->countdown==0)
    {
    taskBuf_group_and_run(buf);
    fprintf(out,"Thread %ul sending signal @ cond %p \n",pthread_self(),buf->finished_cv);
    pthread_cond_broadcast(&buf->finished_cv);
    }
  else
    {
    fprintf(out,"Thread %ul wait @ cond %p \n",pthread_self(),buf->finished_cv);
    // this sends the thread to sleep and releases the mutex. When the signal is
    // received, the mutex is locked again
    pthread_cond_wait(&buf->finished_cv,&buf->lock_mx);
    }
  // release the mutex so the buffer can be used for the next iteration.
  pthread_mutex_unlock(&buf->lock_mx);
  }
  
// once the buffer is full, group the tasks into blocks and enqueue them.
// This function must not be called by multiple threads at a time, which 
// is why we encapsulate it in an openMP lock in taskBuf_exec above.
void taskBuf_group_and_run(taskBuf_t* buf)
  {
  int i,pos1,pos2;
  
  fprintf(out,"control thread %lu starting jobs\n",pthread_self());
  fprintf(out,"job types: ");
  for (i=0;i<buf->njobs;i++)
    {
    fprintf(out," %d",buf->jobType[i]);
    }
  fprintf(out,"\n");
  
  for (i=0;i<njobTypes;i++)
    {
    buf->ghost_args[i]->n=0;
    }
  for (i=0;i<buf->njobs;i++)
      {
      pos1=buf->jobType[i]; // which op is requested? (rndX, incX or divX)
      pos2=buf->ghost_args[pos1]->n; // how many of this op have been requested already?
      buf->ghost_args[pos1]->val[pos2]=buf->jobArg[i];
      buf->ghost_args[pos1]->n++;
      }
  for (i=0;i<njobTypes;i++)
    {
    if (buf->ghost_args[i]->n>0)
      {
      fprintf(out,"enqueue job type %d on %d workers for %d values\n",
        i,nworkers,buf->ghost_args[i]->n);
      ghost_task_add(buf->ghost_task[i]);
      }
    }
  buf->countdown=buf->njobs;
  //TODO - this is a waitall, we should actually just wait for the task that solved
  //       our specific request (e.g. only for incX etc.)
  for (i=0;i<njobTypes;i++)
    {
    if (buf->ghost_args[i]->n>0)
      {
      fprintf(out,"wait for job type %d\n",i);
      ghost_task_wait(buf->ghost_task[i]);
      }
    }  
  return;
  }

////////////////////////////////////////////////////////////////////////////////////////////
// implementation of the algorithm                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
int n;
int col;
int *v;
taskBuf_t *taskBuf;
} mainArg_t;

// this is a long-running job (an "algorithm") executed by a control thread.
// Inside, work that has to be done is put in a task buffer and than bundled
// between control threads and put into the ghost queue.
void* fill_vector(void* v_mainArg)
  {
  int i,n,col;
  int* v;
  taskBuf_t* taskBuf;
  mainArg_t *mainArg = (mainArg_t*)v_mainArg;

  n = mainArg->n;
  col = mainArg->col;
  v = mainArg->v;
  taskBuf=mainArg->taskBuf;
  
  // these are all blocking function calls right now
  v[0]=0;// next one will be randomized
  for (i=1; i<n;i++)
    {
    v[i]=v[i-1];
    if (v[i-1]==0 || v[i-1]==1)
      {
      taskBuf_exec(taskBuf,&v[i],col,JOBTYPE_RNDX);
      }
    else if (v[i-1]%2==0)
      {
      taskBuf_exec(taskBuf,&v[i],col,JOBTYPE_DIVX);
      }
    else
      {
      taskBuf_exec(taskBuf,&v[i],col,JOBTYPE_INCX);
      }
    }
  return;
  }

  
int main(int argc, char** argv)
  {
  int rank, num_proc;
  int ierr;
  int i,j;
  taskBuf_t *taskBuf;
  mainArg_t mainArg[ntasks];
#ifdef DONT_PIN_CONTROL_THREADS
  pthread_t controlThread[ntasks];
#else
  ghost_task_t *controlTask[ntasks];
#endif
  
  int* vectors = (int*)malloc(ntasks*ndim*sizeof(int));
  
  comm_ptr_t comm_world;
  
  out = stderr;
  
  PHIST_CHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);
  // avoid duplicate init call
#ifndef PHIST_KERNEL_LIB_GHOST
  // initialize ghost queue
  ghost_init(argc,argv);
#endif
  PHIST_CHK_IERR(phist_comm_create(&comm_world,&ierr),ierr);
 
#ifdef DONT_PIN_CONTROL_THREADS

nworkers=ghost_thpool->nThreads;

#else
// we would actually like to run the control threads without pinning
// and allow having all cores for doing the work, but if we do that 
// in the current ghost implementation, the jobs in the ghost queue 
// never get executed because ntasks cores are reserved for the     
// control threads.
nworkers=ghost_thpool->nThreads - ntasks;

#endif

PHIST_CHK_IERR(taskBuf_create(&taskBuf,ntasks,nworkers,&ierr),ierr);


for (i=0;i<ntasks;i++)
  {
  mainArg[i].n=ndim;
  mainArg[i].col=i;
  mainArg[i].v=vectors+i*ndim;
  mainArg[i].taskBuf=taskBuf;
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
for (i=0;i<ntasks;i++)
  {
  pthread_join(controlThread[i],NULL);
  }
#endif
// single thread prints the results
for (i=0;i<ndim;i++)
  {
  for (j=0;j<ntasks;j++)
    {
    fprintf(out,"\t%d",vectors[j*ndim+i]);
    }
  fprintf(out,"\n");
  }
  
  PHIST_CHK_IERR(phist_kernels_finalize(&ierr),ierr);
  // avoid duplicate finish call
#ifndef PHIST_KERNEL_LIB_GHOST
  // finalize ghost queue
  ghost_finish();
#endif  
  }
