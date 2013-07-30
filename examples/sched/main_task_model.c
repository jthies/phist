#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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

////////////////////////////////////////////////////////////////////////////////////////
// algorithm-specific functionality                                                   //
////////////////////////////////////////////////////////////////////////////////////////

// argument for rndX, incX and divX functions operating on an array of int pointers.
typedef struct
  {
  int n;
  int **val;
  } argList_t;

void argList_create(argList_t** args)
  {
  *args = (argList_t*)malloc(sizeof(argList_t));
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

#pragma omp parallel
  {
#pragma omp single
  fprintf(stdout,"rndX executed by %d threads on %d values\n",omp_get_num_threads(),arg->n);  
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

#pragma omp parallel
  {
#pragma omp single
  fprintf(stdout,"incX executed by %d threads on %d values\n",
        omp_get_num_threads(),arg->n);  
#pragma omp for
  for (i=0;i<arg->n;i++)
    *(arg->val[i])++;
  }
  return NULL;
  }


// divide all entries in an array by two
static void *divX(void* v_arg)
  {
  int i;
  argList_t* arg = (argList_t*)v_arg;

#pragma omp parallel
  {
#pragma omp single
  fprintf(stdout,"divX executed by %d threads on %d values\n", 
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
omp_lock_t lock; // for controlling access to the shared object.

int *jobType;
int **jobArg;

ghost_task_t **ghost_task;
argList_t **ghost_args;

} taskBuf_t;

void taskBuf_create(taskBuf_t** buf, int num_tasks, int* ierr);
void taskBuf_exec(taskBuf_t* buf, int* arg, int task_id, int taskFlag);
void taskBuf_group_and_enqueue(taskBuf_t* buf);
void taskBuf_wait_all(taskBuf_t* buf);
void taskBuf_delete(taskBuf_t* buf, int* ierr);

// create a new task buffer object
void taskBuf_create(taskBuf_t** buf, int num_tasks, int* ierr)
  {
  int i;
  *ierr=0;
  *buf = (taskBuf_t*)malloc(sizeof(taskBuf_t));
  (*buf)->njobs=num_tasks;
  (*buf)->countdown=num_tasks;
  (*buf)->jobType=(int*)malloc(num_tasks*sizeof(int));
  (*buf)->jobArg=(int**)malloc(num_tasks*sizeof(int*));
  omp_init_lock(&((*buf)->lock));

  (*buf)->ghost_args=(argList_t**)malloc(njobTypes*sizeof(argList_t*));
  (*buf)->ghost_task=(ghost_task_t**)malloc(njobTypes*sizeof(ghost_task_t*));

  for (i=0;i<njobTypes;i++)
    {
    argList_create(&(*buf)->ghost_args[i]);
    }
  
  (*buf)->ghost_task[JOBTYPE_RNDX] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &rndX, 
        (void*)((*buf)->ghost_args[JOBTYPE_RNDX]),GHOST_TASK_DEFAULT);

  (*buf)->ghost_task[JOBTYPE_INCX] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &incX, 
        (void*)((*buf)->ghost_args[JOBTYPE_INCX]),GHOST_TASK_DEFAULT);

  (*buf)->ghost_task[JOBTYPE_DIVX] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &divX, 
        (void*)((*buf)->ghost_args[JOBTYPE_DIVX]),GHOST_TASK_DEFAULT);
  
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
  omp_destroy_lock(&buf->lock);
  for (i=0;i<njobTypes;i++)
    {
    ghost_task_destroy(buf->ghost_task[i]);
    argList_delete(buf->ghost_args[i]);
    }
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
  omp_set_lock(&buf->lock);
  buf->jobType[task_id]=taskFlag;
  buf->jobArg[task_id]=arg;
  buf->countdown--;
  fprintf(stdout,"control thread %d exec job %d on column %d, countdown=%d\n",
        omp_get_thread_num(), taskFlag, task_id, buf->countdown);
  if (buf->countdown==0) 
    {
    taskBuf_group_and_enqueue(buf);
    }
  omp_unset_lock(&buf->lock);
  
  // wait for the grouped tasks to finish (there are typically
  // less tasks enqueued than in the task buffer because of the
  // grouping, for instance: 1x rndX, 2x incX, 1x divX => 3 tasks)
  taskBuf_wait_all(buf);
  // the correct result is now in *arg.
  }
  
// once the buffer is full, group the tasks into blocks and enqueue them.
// This function must not be called by multiple threads at a time, which 
// is why we encapsulate it in an openMP lock in taskBuf_exec above.
void taskBuf_group_and_enqueue(taskBuf_t* buf)
  {
  int i,pos1,pos2;
  
  fprintf(stdout,"control thread %d starting jobs\n",omp_get_thread_num());
  fprintf(stdout,"job types: ");
  for (i=0;i<buf->njobs;i++)
    {
    fprintf(stdout," %d",buf->jobType[i]);
    }
  fprintf(stdout,"\n");
  
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
      ghost_task_add(buf->ghost_task[i]);
      }
    }
  buf->countdown=buf->njobs;
  return;
  }

void taskBuf_wait_all(taskBuf_t* buf)
  {
  // note: ghost_task_waitall alone won't do the trick here because the
  // jobs are typically still in the task buffer and not yet in the
  // ghost queue.
  ghost_task_waitall();
  fprintf(stdout,"thread %d at barrier, waiting for team of %d\n",
        omp_get_thread_num(),omp_get_num_threads());
  }

// this is a long-running job (an "algorithm") executed by a control thread.
// Inside, work that has to be done is put in a task buffer and than bundled
// between control threads and put into the ghost queue.
void fill_column(int* v0, int n, int col, taskBuf_t* taskBuf)
  {
  int i;
  int* v = v0 + n*col;
  fprintf(stdout,"control thread %d fills vector %d\n",
        omp_get_thread_num(), col);
  // these are all blocking function calls right now
  v[0]=0;// next one will be randomized
  for (i=1; i<n;i++)
    {
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

  int* vectors = (int*)malloc(ntasks*ndim*sizeof(int));
  
  comm_ptr_t comm_world;
  
  PHIST_CHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);
  // avoid duplicate init call
#ifndef PHIST_KERNEL_LIB_GHOST
  // initialize ghost queue
  ghost_init(argc,argv);
#endif
  PHIST_CHK_IERR(phist_comm_create(&comm_world,&ierr),ierr);
 
PHIST_CHK_IERR(taskBuf_create(&taskBuf,ntasks,&ierr),ierr);

// create the control threads. A single thread creates ntasks tasks, 
// which are executed by anyone. Each task puts small jobs in a buffer,
// which are bundled and enqueued in the ghost queue. Nested parallelism
// is allowed in OpenMP 3.0, so whenever tasks are executed by ghost, 
// a single control thread creates a team of worker threads to do the work.
omp_set_num_threads(ntasks);
#pragma omp parallel
{
fprintf(stdout,"this is control thread %d of %d\n",
        omp_get_thread_num(),omp_get_num_threads());
#pragma omp single
  {
  fprintf(stdout,"thread %d of %d starts all tasks.\n",
          omp_get_thread_num(),omp_get_num_threads());    
#pragma omp task untied
    fill_column(vectors,ndim,0, taskBuf);
#pragma omp task untied
    fill_column(vectors,ndim,1, taskBuf);
#pragma omp task untied
    fill_column(vectors,ndim,2, taskBuf);
#pragma omp task untied
    fill_column(vectors,ndim,3, taskBuf);
    //TODO - the loop variant doesn't work, multiple threads grab a single column
    //       (the i variable is not correct in the private data scope created?)
//  for (i=0;i<ntasks;i++)
//    {
//#pragma omp task untied private(i)
//    fill_column(vectors,ndim,i, taskBuf);
//    }
  }
#pragma omp taskwait
}
  
// single thread prints the results
for (i=0;i<ndim;i++)
  {
  for (j=0;j<ntasks;j++)
    {
    fprintf(stdout,"\t%d",vectors[j*ndim+i]);
    }
  fprintf(stdout,"\n");
  }
  
  PHIST_CHK_IERR(phist_kernels_finalize(&ierr),ierr);
  // avoid duplicate finish call
#ifndef PHIST_KERNEL_LIB_GHOST
  // finalize ghost queue
  ghost_finish();
#endif  
  }
