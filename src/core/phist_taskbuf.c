#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "phist_taskbuf.h"
#include "phist_macros.h"

void argList_create(argList_t** args, int ntasks, int* ierr)
  {
  *ierr=0;
  *args = (argList_t*)malloc(sizeof(argList_t));
  (*args)->nthreads=1;
  (*args)->n=0;
  (*args)->arg=(void**)malloc(ntasks*sizeof(void*));
  }

void argList_delete(argList_t* args, int* ierr)
  {
  *ierr=0;
  free(args->arg);
  free(args);
  }


// create a new task buffer object
void taskBuf_create(taskBuf_t** buf, int num_controllers, int* ierr)
  {
  int i;
  *ierr=0;
  *buf = (taskBuf_t*)malloc(sizeof(taskBuf_t));
  
  (*buf)->ncontrollers=num_controllers;
  (*buf)->countdown=num_controllers;
  (*buf)->opType=(int*)malloc(num_controllers*sizeof(int));
  (*buf)->opArg=(void**)malloc(num_controllers*sizeof(int*));

  pthread_mutex_init(&((*buf)->lock_mx),NULL);
  pthread_cond_init(&(*buf)->finished_cv,NULL);
  
  (*buf)->nops = 0;
  (*buf)->nops_max = 8;

  (*buf)->ghost_args=(argList_t**)malloc((*buf)->nops_max*sizeof(argList_t*));
  (*buf)->ghost_task=(ghost_task_t**)malloc((*buf)->nops_max*sizeof(ghost_task_t*));

  return;
  }

// create a new task buffer object
void taskBuf_delete(taskBuf_t* buf, int* ierr)
  {
  int i;
  *ierr=0;
  //TODO - maybe we should force execution/wait here
  free(buf->opType);
  free(buf->opArg);
  for (i=0;i<buf->nops;i++)
    {
    ghost_task_destroy(buf->ghost_task[i]);
    PHIST_CHK_IERR(argList_delete(buf->ghost_args[i],ierr),*ierr);
    }
  pthread_mutex_destroy(&buf->lock_mx);
  pthread_cond_destroy(&buf->finished_cv);
  free(buf->ghost_task);
  free(buf->ghost_args);
  free(buf);
  return;
  }

//! add a function to the buffer so that taskBuf_add can be used to request
//! that operation. The key is generated automatically by this function and can be used
//! to request the operation later on.
void taskBuf_add_op(taskBuf_t* buf, void* (*fcn)(argList_t*),
        int num_workers, int* op_key, int* ierr)
  {
  *ierr=0;
  if (buf->nops==buf->nops_max)
    {
    PHIST_OUT(1,"%d ops in task buffer, re-allocating memory",buf->nops+1);
    buf->nops_max += 8;
    buf->ghost_args=(argList_t**)realloc(buf->ghost_args,buf->nops_max*sizeof(argList_t*));
    buf->ghost_task=(ghost_task_t**)realloc(buf->ghost_task,buf->nops_max*sizeof(ghost_task_t*));
    }
  *op_key = buf->nops;
  buf->nops++;
  PHIST_CHK_IERR(argList_create(&buf->ghost_args[*op_key], buf->ncontrollers,ierr),*ierr);
  buf->ghost_args[*op_key]->nthreads = num_workers;
  buf->ghost_args[*op_key]->nthreads = num_workers;
  //TODO - properly state in which queue (LD) the job should be put
  if (num_workers==ghost_thpool->nThreads)
    {
    buf->ghost_task[*op_key] = ghost_task_init(GHOST_TASK_FILL_ALL, 0, (void*(*)(void*))fcn, 
        (void*)buf->ghost_args[*op_key],GHOST_TASK_DEFAULT);
    }
  else
    {
    buf->ghost_task[*op_key] = ghost_task_init(num_workers, 0, (void*(*)(void*))fcn, 
        (void*)(buf->ghost_args[*op_key]),GHOST_TASK_DEFAULT);
    }
  }

// a job is put into the task buffer and executed once the
// buffer is full. This function returns only when the task 
// is finished.
void taskBuf_add(taskBuf_t* buf, void* arg, int task_id, int op_key,int* ierr)
  {
  *ierr=0;
  // lock the buffer while putting in jobs and executing them
  pthread_mutex_lock(&buf->lock_mx);
  buf->opType[task_id]=op_key;
  buf->opArg[task_id]=arg;
  buf->countdown--;
  PHIST_OUT(1,"control thread %lu request job %d on column %d, countdown=%d\n",
        pthread_self(), op_key, task_id, buf->countdown);
  if (buf->countdown==0)
    {
    taskBuf_flush(buf);
    PHIST_OUT(1,"Thread %lu sending signal @ cond %p \n",pthread_self(),buf->finished_cv);
    pthread_cond_broadcast(&buf->finished_cv);
    }
  else
    {
    PHIST_OUT(1,"Thread %lu wait @ cond %p \n",pthread_self(),buf->finished_cv);
    // this sends the thread to sleep and releases the mutex. When the signal is
    // received, the mutex is locked again
    pthread_cond_wait(&buf->finished_cv,&buf->lock_mx);
    }
  // the bundled tasks have been put in the ghost queue.
  // Release the mutex so the buffer can be used for the next iteration.
  // The control thread is responsible for waiting for the task to finish
  // before putting in a new one, otherwise a race condition is created.
  pthread_mutex_unlock(&buf->lock_mx);
  }
  
// once the buffer is full, group the tasks into blocks and enqueue them.
// This function must not be called by multiple threads at a time, which 
// is why we encapsulate it in an openMP lock in taskBuf_add above.
void taskBuf_flush(taskBuf_t* buf)
  {
  int i,pos1,pos2;
  
  PHIST_OUT(1,"control thread %lu starting jobs\n",pthread_self());
  PHIST_OUT(1,"job types: ");
  for (i=0;i<buf->ncontrollers;i++)
    {
    PHIST_OUT(1," %d",buf->opType[i]);
    }
  PHIST_OUT(1,"\n");
  
  for (i=0;i<buf->nops;i++)
    {
    buf->ghost_args[i]->n=0;
    }
  for (i=0;i<buf->ncontrollers;i++)
      {
      pos1=buf->opType[i]; // which op is requested? (rndX, incX or divX)
      if (pos1>=0) // otherwise: -1=NO-OP
        {
        pos2=buf->ghost_args[pos1]->n; // how many of this op have been requested already?
        buf->ghost_args[pos1]->arg[pos2]=buf->opArg[i];
        buf->ghost_args[pos1]->n++;
        }
      }
  for (i=0;i<buf->nops;i++)
    {
    if (buf->ghost_args[i]->n>0)
      {
      PHIST_OUT(1,"enqueue job type %d on %d workers for %d values\n",
        i,buf->ghost_args[i]->nthreads,buf->ghost_args[i]->n);
      ghost_task_add(buf->ghost_task[i]);
      }
    }
  buf->countdown=buf->ncontrollers;
#ifdef SYNC_WAIT
  //TODO
  // this should be removed, waiting should be left to 
  // each of the control threads using taskBuf_wait().
  for (i=0;i<buf->nops;i++)
  {
  if (buf->ghost_args[i]->n>0)
    {
    PHIST_OUT(1,"Thread %lu waiting for job type %d\n",pthread_self(),i);
      ghost_task_wait(buf->ghost_task[i]);
    }
  }
#endif
  return;
  }

// if a thread finishes with whatever it is the task buffer has been created    
// for and others should continue using the buffer, the thread can renounce     
// the buffer and others will no longer wait for this thread to add its         
// requests.
void taskBuf_renounce(taskBuf_t* buf, int task_id, int* ierr)
  {
  *ierr=0;
  // lock the buffer while adjusting it to deal without me in the future.
  pthread_mutex_lock(&buf->lock_mx);
  buf->opType[task_id]=PHIST_NOOP;
  buf->opArg[task_id]=NULL;
  buf->ncontrollers--;
  buf->countdown--;
  PHIST_OUT(1,"control thread %lu (task_id %d) leaving the controller team\n",
        pthread_self(), task_id);
  // final flush if necessary
  if (buf->countdown==0)
    {
    taskBuf_flush(buf);
    PHIST_OUT(1,"Thread %lu sending signal @ cond %p \n",pthread_self(),buf->finished_cv);
    pthread_cond_broadcast(&buf->finished_cv);
    }
  pthread_mutex_unlock(&buf->lock_mx);  
  return;
  }

void taskBuf_wait(taskBuf_t* buf, int task_id, int* ierr)
  {
  int jt = buf->opType[task_id];
  *ierr=0;
  if (buf->ghost_args[jt]->n>0)
    {
    PHIST_OUT(1,"Thread %lu waiting for job type %d\n",pthread_self(),jt);
      ghost_task_wait(buf->ghost_task[jt]);
    PHIST_OUT(1,"Thread %lu, job type %d finished.\n",pthread_self(),jt);
    }
  return;
  }

