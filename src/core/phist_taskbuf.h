#ifndef PHIST_TASKBUF_H
#define PHIST_TASKBUF_H

#include <pthread.h>

#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#include "ghost_taskq.h"
#else
#error "ghost is required here."
#endif

#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////////////////////////////////
// task buffer prototype                                                              //
////////////////////////////////////////////////////////////////////////////////////////

//! argument for functions operating on multiple pointers simultaneously
//! (the kind of functions we want to bundle in the task buffer)
typedef struct
  {
  int nthreads;
  int n;
  void **arg;
  } argList_t;

void argList_create(argList_t** args, int ntasks, int* ierr);
void argList_delete(argList_t* args, int* ierr);

// a buffer object, when each control thread has set its opType flag,
// the operations are bundled and enqueued in the ghost queue. The
// members of this struct should never be accessed by the user but
// only by taskBuf_* functions
typedef struct {

int ncontrollers; // total number of threads using the buffer
int nops; // number of different operations the buffer can perform
int nops_max; // number of ops that can be added before memory has to be re-allocated
int countdown; // when this reaches 0, the jobs are launched
pthread_mutex_t lock_mx; // for controlling access to the shared object and the condition 
                         // variable finished_cv.
pthread_cond_t finished_cv; // for waiting for the buffer to be emptied

int *opType; // dimension ncontrollers
void **opArg; // dimension ncontrollers
int *opWorkers; // dimension ncontrollers, number of work threads to be employed for the op

ghost_task_t **ghost_task; // dimension nops
argList_t **ghost_args;    // dimension nops

} taskBuf_t;

//! create a new task buffer for handling certain basic operations in a blocked way
void taskBuf_create(taskBuf_t** buf, int num_controllers, int* ierr);

//! delete the task buffer
void taskBuf_delete(taskBuf_t* buf, int* ierr);

//! add a function to the buffer so that taskBuf_add can be used to request
//! that operation. The key is generated automatically by this function and can be used
//! to request the operation later on.
void taskBuf_add_op(taskBuf_t* buf, void* (*fcn)(argList_t*), 
        int num_workers, int* op_key, int* ierr);

//! put a job request into the task buffer and wait until it is enqueued
void taskBuf_add(taskBuf_t* buf, void* arg, int task_id, int op_key, int* ierr);

//! flush the task buffer, e.g. group the tasks together and put them in
//! the ghost queue, then signal any waiting threads.
void taskBuf_flush(taskBuf_t* buf);

//! wait for the task that you put in the buffer to be finished
void taskBuf_wait(taskBuf_t* buf, int task_id, int* ierr);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
