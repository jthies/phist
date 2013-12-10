#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

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
  (*args)->shared_arg=NULL;
  (*args)->in_arg=(const void**)malloc(ntasks*sizeof(void*));
  (*args)->out_arg=(void**)malloc(ntasks*sizeof(void*));
  (*args)->id=(int*)malloc(ntasks*sizeof(int));
  }

void argList_delete(argList_t* args, int* ierr)
  {
  *ierr=0;
  free(args->in_arg);
  free(args->out_arg);
  free(args->id);
  free(args);
  }


// create a new task buffer object
void taskBuf_create(taskBuf_t** buf, int num_controllers, int* ierr)
  {
  int i;
  *ierr=0;
  *buf = (taskBuf_t*)malloc(sizeof(taskBuf_t));
  
  (*buf)->ncontrollers=num_controllers;
  (*buf)->nactive=num_controllers;
  (*buf)->countdown=num_controllers;
  (*buf)->opType=(int*)malloc(num_controllers*sizeof(int));
  (*buf)->waitFor=(int*)malloc(num_controllers*sizeof(int));
  (*buf)->array_mx=(pthread_mutex_t*)malloc(num_controllers*sizeof(pthread_mutex_t));
  for (i=0;i<num_controllers;i++)
    {
    (*buf)->opType[i]=PHIST_OP_NULL;
    (*buf)->waitFor[i]=PHIST_OP_NULL;
    pthread_mutex_init(&((*buf)->array_mx[i]),NULL);
    }
  (*buf)->op_inArgs=(const void**)malloc(num_controllers*sizeof(void*));
  (*buf)->op_outArgs=(void**)malloc(num_controllers*sizeof(void*));

  pthread_mutex_init(&((*buf)->lock_mx),NULL);
  pthread_cond_init(&(*buf)->flushed_cv,NULL);
  
  (*buf)->nops = 0;
  (*buf)->nops_max = 32; // that should be enough for all practical purposes, but
                         // if it is ever exceeded we reallocate in add_op

  (*buf)->ghost_args=(argList_t**)malloc(((*buf)->nops_max+1)*sizeof(argList_t*));
  (*buf)->ghost_task=(ghost_task_t**)malloc(((*buf)->nops_max+1)*sizeof(ghost_task_t*));

  return;
  }

// create a new task buffer object
void taskBuf_delete(taskBuf_t* buf, int* ierr)
  {
  int i;
  *ierr=0;
  //TODO - maybe we should force execution/wait here
  free(buf->opType);
  free(buf->waitFor);
  free(buf->op_inArgs);
  free(buf->op_outArgs);
  for (i=1;i<=buf->nops;i++)
    {
    ghost_task_destroy(buf->ghost_task[i]);
    PHIST_CHK_IERR(argList_delete(buf->ghost_args[i],ierr),*ierr);
    pthread_mutex_destroy(&buf->array_mx[i]);
    }
  pthread_mutex_destroy(&buf->lock_mx);
  pthread_cond_destroy(&buf->flushed_cv);
  free(buf->ghost_task);
  free(buf->ghost_args);
  free(buf);
  return;
  }

//! add a function to the buffer so that taskBuf_add can be used to request
//! that operation. The key is generated automatically by this function and can be used
//! to request the operation later on.
void taskBuf_add_op(taskBuf_t* buf, void* (*fcn)(argList_t*), void* shared_arg,
        int num_workers, int* op_key, int* ierr)
  {
  *ierr=0;
  if (buf->nops==buf->nops_max)
    {
    PHIST_OUT(2,"%d ops in task buffer, re-allocating memory",buf->nops+1);
    buf->nops_max *= 2;
    buf->ghost_args=(argList_t**)realloc(buf->ghost_args,(buf->nops_max+1)*sizeof(argList_t*));
    buf->ghost_task=(ghost_task_t**)realloc(buf->ghost_task,(buf->nops_max+1)*sizeof(ghost_task_t*));
    }
  // note that slot 0 always remains empty, this is to be consistent with PHIST_OP_NULL=0
  buf->nops++;
  *op_key = buf->nops;
  PHIST_CHK_IERR(argList_create(&buf->ghost_args[*op_key], buf->ncontrollers,ierr),*ierr);
  buf->ghost_args[*op_key]->nthreads = num_workers;
  buf->ghost_args[*op_key]->shared_arg=shared_arg;
  buf->ghost_args[*op_key]->nthreads = num_workers;
  //TODO - properly state in which queue (LD) the job should be put
  if (num_workers==ghost_thpool->nThreads)
    {
    buf->ghost_task[*op_key] = ghost_task_init(GHOST_TASK_FILL_ALL, GHOST_TASK_LD_ANY, 
        (void*(*)(void*))fcn, (void*)buf->ghost_args[*op_key],GHOST_TASK_DEFAULT);
    }
  else
    {
    // TODO - correct usage of locality domains
    buf->ghost_task[*op_key] = ghost_task_init(num_workers, GHOST_TASK_LD_ANY, 
        (void*(*)(void*))fcn, (void*)(buf->ghost_args[*op_key]),GHOST_TASK_DEFAULT);
    }
  }

// a job is put into the task buffer and enqueued once the
// buffer is full. This function returns immediately. Before
// calling it again, you have to call taskBuf_wait to make  
// sure the task has been enqueued and is finished.
void taskBuf_add(taskBuf_t* buf, const void* in_arg, 
                                       void* out_arg, 
                             int task_id, int op_key, int* ierr)
  {
  *ierr=0;
  // first make sure that the previous operation has been waited for
  // We need to actually wait for completion of the ghost task here.
  // Otherwise we could get that a ghost task is still in the
  // queue or running and the ghost_task_t is overwritten by  
  // the next taskBuf_flush. On the other hand, the buffer is only  
  // flushed if all participating threads have put in something,    
  // which by this mechanism makes sure that all ghost tasks have   
  // finished before the next round is started.
  pthread_mutex_lock(&buf->array_mx[task_id]);
  if (buf->opType[task_id]!=PHIST_OP_NULL)
    {
    PHIST_OUT(2,"forcing thread %lu to wait for previous flush\n",pthread_self());
    pthread_mutex_unlock(&buf->array_mx[task_id]);
    PHIST_CHK_IERR(taskBuf_wait(buf,task_id,ierr),*ierr);
    pthread_mutex_lock(&buf->array_mx[task_id]);
    }
  if (buf->waitFor[task_id]!=PHIST_OP_NULL)
    {
    PHIST_OUT(2,"forcing thread %lu to wait for previous job to finish\n",pthread_self());
    pthread_mutex_unlock(&buf->array_mx[task_id]);
    PHIST_CHK_IERR(taskBuf_wait(buf,task_id,ierr),*ierr);
    pthread_mutex_lock(&buf->array_mx[task_id]);
    }

  // lock the buffer while putting in jobs and executing them
  buf->opType[task_id]=op_key;
  buf->waitFor[task_id]=PHIST_OP_NULL; // indicate that we need a flush before ghost_task_wait
  buf->op_inArgs[task_id]=in_arg;
  buf->op_outArgs[task_id]=out_arg;
  pthread_mutex_unlock(&buf->array_mx[task_id]);
  
  pthread_mutex_lock(&buf->lock_mx);
  buf->countdown--;
  PHIST_OUT(2,"control thread %lu request job %d on column %d, countdown=%d\n",
        pthread_self(), op_key, task_id, buf->countdown);
  // the first thread who observes that the buffer is full
  // flushes it, which means that the operations are bundled
  // and put in the ghost queue. From this moment the buffer
  // can be used for the next iteration. It is the responsibility
  // of the control threads to do a taskBuf_wait before reusing the
  // buffer.
  if (buf->countdown==0)
    {
    taskBuf_flush(buf);
    }
  pthread_mutex_unlock(&buf->lock_mx);
  }
  
// once the buffer is full, group the tasks into blocks and enqueue them.
// This function must not be called by multiple threads at a time, which 
// is why we encapsulate it in a pthread lock in taskBuf_add above. This 
// function is NOT intended for use by anyone (consider it private)
void taskBuf_flush(taskBuf_t* buf)
  {
  int i,pos1,pos2;
  
  PHIST_OUT(2,"control thread %lu starting jobs\n",pthread_self());
  PHIST_OUT(2,"job types: ");
  
  // skip op 0, which is PHIST_OP_NULL ("do nothing")
  for (i=1;i<=buf->nops;i++)
    {
    buf->ghost_args[i]->n=0;
    }
  for (i=0;i<buf->ncontrollers;i++)
    {
    pthread_mutex_lock(&buf->array_mx[i]);
    PHIST_OUT(2," %d",buf->opType[i]);
    PHIST_OUT(2,"\n");
    pos1=buf->opType[i]; // which op is requested? (rndX, incX or divX)
    buf->opType[i]=PHIST_OP_NULL; // next op can be requested (buffer has been flushed)
    buf->waitFor[i]=pos1;
    pthread_mutex_unlock(&buf->array_mx[i]);
    if (pos1>0) // otherwise: 0=PHIST_OP_NULL
      {
      pos2=buf->ghost_args[pos1]->n; // how many of this op have been requested already?
      buf->ghost_args[pos1]->id[pos2]=i; // some operators find it interesting to know
                                         // which column they operate on.
                                         // (TODO: do we actually need this?)
      buf->ghost_args[pos1]->in_arg[pos2]=buf->op_inArgs[i];
      buf->ghost_args[pos1]->out_arg[pos2]=buf->op_outArgs[i];
      buf->ghost_args[pos1]->n++;
      }
    }
  for (i=1;i<=buf->nops;i++)
    {
    PHIST_OUT(3,"job type %d: %d requests\n",i,buf->ghost_args[i]->n);
    if (buf->ghost_args[i]->n>0)
      {
      PHIST_OUT(2,"enqueue job type %d on %d workers for %d values\n",
        i,buf->ghost_args[i]->nthreads,buf->ghost_args[i]->n);
      ghost_task_add(buf->ghost_task[i]);
      }
    }
  buf->countdown=buf->nactive;
  PHIST_OUT(2,"Thread %lu sending signal @ cond %p \n",pthread_self(),&buf->flushed_cv);
  pthread_cond_broadcast(&buf->flushed_cv);
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
  pthread_mutex_lock(&buf->array_mx[task_id]);
  buf->opType[task_id]=PHIST_OP_NULL;
  buf->op_inArgs[task_id]=NULL;
  buf->op_outArgs[task_id]=NULL;
  pthread_mutex_unlock(&buf->array_mx[task_id]);
  pthread_mutex_lock(&buf->lock_mx);
  buf->nactive--;
  buf->countdown--;
  PHIST_OUT(2,"control thread %lu (task_id %d) leaving the controller team\n",
        pthread_self(), task_id);
  // final flush if necessary
  if (buf->countdown==0)
    {
    taskBuf_flush(buf);
    }
  pthread_mutex_unlock(&buf->lock_mx);  
  return;
  }

void taskBuf_wait(taskBuf_t* buf, int task_id, int* ierr)
  {
  int jt;
  int wf;
  *ierr=0;
  // lock the array position so we don't miss the taskBuf_flush
  pthread_mutex_lock(&buf->array_mx[task_id]);
  jt = buf->opType[task_id];
  wf = buf->waitFor[task_id];

  // no operation submitted -> return immediately
  if (jt==PHIST_OP_NULL && wf==PHIST_OP_NULL) 
    {
    pthread_mutex_unlock(&buf->array_mx[task_id]);
    return; 
    }


  // wait until the buffer has been flushed
  if (jt!=PHIST_OP_NULL)
    {
    pthread_mutex_lock(&buf->lock_mx);
    PHIST_OUT(2,"Thread %lu waiting for flush @ %p\n",pthread_self(),&buf->flushed_cv);
    // if we don't release the array position lock, taskBuf_flush will
    // create a deadlock
    pthread_mutex_unlock(&buf->array_mx[task_id]);
    pthread_cond_wait(&buf->flushed_cv,&buf->lock_mx);
    pthread_mutex_lock(&buf->array_mx[task_id]);
    pthread_mutex_unlock(&buf->lock_mx);
    }
  // by now we are sure that someone has called taskBuf_flush and
  // the requested operation is in the queue.
  
  // it is safe to keep the array position task_id locked until the
  // job has been executed because no flush can occur while we're  
  // still here (that is, until this thread announces it's next op)
  wf = buf->waitFor[task_id];

  if (wf!=PHIST_OP_NULL)
    {
    PHIST_OUT(2,"Thread %lu waiting for job type %d\n",pthread_self(),wf);
    ghost_task_wait(buf->ghost_task[wf]);
    buf->waitFor[task_id]=PHIST_OP_NULL;
    PHIST_OUT(2,"Thread %lu, job type %d finished.\n",pthread_self(),jt);
    }
  pthread_mutex_unlock(&buf->array_mx[task_id]);
  return;
  }

