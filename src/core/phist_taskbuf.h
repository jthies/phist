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

#define PHIST_OP_NULL 0

//! implementation of a task buffer that bundles operations which    
//! are memory intensive because of a shared operation component.    
//! The main applications in phist are the shifted operator (A-sjI)X 
//! and the projected preconditioning M\(I-VV')X, where the columns  
//! in X are managed by different threads.                           

//! For an example of how the task buffer is used, see 
//! examples/sched/main_task_model.c

//! argument for functions operating on multiple pointers simultaneously
//! (the kind of functions we want to bundle in the task buffer).       
//!                                                                     
//! The canonical example is compute Y[i] = A*X[i] for n (multi-)vectors
//! X, Y and a single operator A, using up to nthreads worker threads.  
//! Here A would be the shared_arg, X[i] the in_args and Y[i] the       
//! out_args.
typedef struct
  {
  //! \name common arguments useful for any function used with a buffer
  //!@{
  int nthreads; //! how many (OpenMP) worker threads are reserved for this task?
  int n; //! how many arguments (vectors, columns or whatever) should you work on?
  int* id; //! stores the ID's (for instance column index) TODO - remove this, I think
  //@}
  const void* shared_arg; //! an input argument required for all vectors
  const void* *in_arg; //! array of input arguments passed to the function (for instance vectors)
  void** out_arg; //! array of output arguments passed to the function
  int ierr;
  } argList_t;

void argList_create(argList_t** args, int ntasks, int* ierr);
void argList_delete(argList_t* args, int* ierr);

// a buffer object, when each control thread has set its opType flag,
// the operations are bundled and enqueued in the ghost queue. The
// members of this struct should never be accessed by the user but
// only by taskBuf_* functions
typedef struct {

int ncontrollers; // total number of threads using the buffer
int nactive; // number of active controllers (some may renounce the buffer at a point)
int nops; // number of different operations the buffer can perform
int nops_max; // number of ops that can be added before memory has to be re-allocated
int countdown; // when this reaches 0, the jobs are launched
pthread_mutex_t lock_mx; // for controlling access to the shared object and the condition 
                         // variable flushed_cv.
pthread_cond_t flushed_cv; // for waiting for the buffer to be emptied
pthread_mutex_t *array_mx; // for controlling access to the arrays opType and waitFor,
                          // where thread i may access element i, but also others if 
                          // he is the one flushing the buffer.

int *opType; // dimension ncontrollers: which operation is requested by thread i next?
int *waitFor; // dimension ncontrollers: after flush, which op should task i wait for?

// all others are only accessed atomically in taskBuf_flush, so they are protected by lock_mx
const void **op_inArgs; // when a task is added for a certain column, these are the input args
void **op_outArgs;      // and these the output args provided by the thread. They may have to be
                         // reordered before put in the ghost queue.

ghost_task_t **ghost_task; // dimension nops: task to bundle op i and enqueue it in ghost
argList_t **ghost_args;    // dimension nops: array of args passed into task i

} taskBuf_t;

//! create a new task buffer for handling certain basic operations in a blocked way
void taskBuf_create(taskBuf_t** buf, int num_controllers, int* ierr);

//! delete the task buffer
void taskBuf_delete(taskBuf_t* buf, int* ierr);

//! add a function to the buffer so that taskBuf_add can be used to request
//! that operation. The key is generated automatically by this function and can be used
//! to request the operation later on.
void taskBuf_add_op(taskBuf_t* buf, void* (*fcn)(argList_t*), 
        void* shared_arg, int num_workers, int* op_key, int* ierr);

//! put a job request into the task buffer. This is a non-blocking call that
//! returns immediately. If taskBuf_wait has not been called since the last 
//! call to taskBuf_add with this task_id, it is enforced in this function. 
//! Typically the user should call taskBuf_wait before another taskBuf_add  
//! with the same task_id. The special argument op_key=PHIST_OP_NULL can be 
//! given to this function to indicate that for this task_id no operation   
//! should be done in this round. If no more operations are desired for this
//! task_id, taskBuf_renounce should be called instead of taskBuf_add.
void taskBuf_add(taskBuf_t* buf, const void* in_arg, void* out_arg, 
                             int task_id, int op_key, int* ierr);

//! flush the task buffer, i.e. group the tasks together and put them in
//! the ghost queue, then signal any waiting threads. This is a private 
//! function that should never be called from the outside. It's call
//! must be protected by buf->lock_mx.
void taskBuf_flush(taskBuf_t* buf);

//! if a thread finishes with whatever it is the task buffer has been created    
//! for and others should continue using the buffer, the thread can renounce     
//! the buffer and others will no longer wait for this thread to add its         
//! requests.
void taskBuf_renounce(taskBuf_t* buf, int task_id, int* ierr);


//! wait for the task that you put in the buffer to be finished.
//! This includes (a) waiting for the buffer to be flushed, and 
//! (b) waiting for the job to be run and finished via the ghost
//! queue. This is a synchronization point between all threads  
//! active in the buffer.
void taskBuf_wait(taskBuf_t* buf, int task_id, int* ierr);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
