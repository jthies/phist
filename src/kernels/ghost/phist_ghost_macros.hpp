#ifndef PHIST_GHOST_MACROS_HPP
#define PHIST_GHOST_MACROS_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"

#include <ghost.h>
#include <ghost/machine.h>


#ifdef PHIST_HAVE_CXX11_LAMBDAS

// glue between C++11 Lambda-functions and plain old C function pointers
template<typename LFunc>
void phist_execute_lambda_as_ghost_task(ghost_task_t **task, LFunc context, int* iflag, bool async)
{
  // ghost requires a task function of the form "void* f(void*)"
  // a Lambda without capture can be converted automatically to a plain old C function pointer,
  // but the context goes out of scope when this function is left (so we copy it to the heap and destroy it here!)
  auto void_lambda_caller = [] (void* context) -> void* {
    LFunc* pContext = (LFunc*)context;  // get the pointer to the context-lambda-object
    (*pContext)();                      // actually run the target function
    delete pContext;                    // delete the copy of the lambda-function-object
    return NULL;                        // we do not return anything
  };

  // actually create the task and copy the context to the heap
  int nThreads = async ? 0 : GHOST_TASK_FILL_ALL;
  ghost_task_flags_t flags = async ? GHOST_TASK_NOT_PIN : GHOST_TASK_DEFAULT;
  PHIST_CHK_GERR(ghost_task_create(task, nThreads, 0, void_lambda_caller, (void*) new LFunc(context), flags, NULL, 0), *iflag);

  // we do not use the task hierarchy - all tasks are children of the main task!
  // So overwrite the tasks parent!
  ghost_task_t *parentTask = NULL;
  ghost_task_cur(&parentTask);
  if( parentTask )
    while(parentTask->parent)
      parentTask = parentTask->parent;
//  (*task)->parent = parentTask;

PHIST_DEB("enqueuing C++11-lambda as GHOST task and waiting for it\n");
  PHIST_CHK_GERR(ghost_task_enqueue(*task), *iflag);
  if( async )
    return;

  PHIST_CHK_GERR(ghost_task_wait(*task), *iflag);
  ghost_task_destroy(*task);
  *task = NULL;
}

inline void phist_wait_ghost_task(ghost_task_t** task, int* iflag)
{
  PHIST_CHK_IERR(*iflag = (*task != NULL) ? PHIST_SUCCESS : PHIST_WARNING, *iflag);

  PHIST_CHK_GERR(ghost_task_wait(*task), *iflag);
  ghost_task_destroy(*task);
  *task = NULL;
}

// some helpful macros

#define PHIST_TASK_BEGIN(taskName) ghost_task_t* taskName = NULL; phist_execute_lambda_as_ghost_task(&taskName, [&]() -> void {
#define PHIST_TASK_END(task_ierr)          }, task_ierr, false);PHIST_CHK_IERR((void)*(task_ierr),*(task_ierr));
#define PHIST_TASK_END_NOWAIT(task_ierr)   }, task_ierr, true );PHIST_CHK_IERR((void)*(task_ierr),*(task_ierr));
#define PHIST_TASK_WAIT(taskName,task_ierr) PHIST_CHK_IERR(phist_wait_ghost_task(&taskName,task_ierr),*(task_ierr));

#else /* PHIST_HAVE_CXX11_LAMBDAS */

#warning "C++11 not supported, not using GHOST tasking mechanism!"

#define PHIST_TASK_BEGIN(taskName)
#define PHIST_TASK_END(task_ierr)
#define PHIST_TASK_END_NOWAIT(task_ierr)
#define PHIST_TASK_WAIT(taskName,task_ierr)

#endif /* PHIST_HAVE_CXX11_LAMBDAS */

#endif /* PHIST_GHOST_MACROS_HPP */
