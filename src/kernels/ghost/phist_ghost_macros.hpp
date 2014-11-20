#ifndef PHIST_GHOST_MACROS_HPP
#define PHIST_GHOST_MACROS_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <ghost.h>
#include <ghost/machine.h>


#ifdef PHIST_HAVE_CXX11_LAMBDAS

// glue between C++11 Lambda-functions and plain old C function pointers
template<typename LFunc>
void phist_execute_lambda_as_ghost_task(LFunc context)
{
  // a Lambda without capture can be converted automatically to a plain old C function pointer
  auto void_lambda_caller = [] (void* context) -> void* {(*(LFunc*)context)(); return NULL;};

  ghost_task_t *task = NULL;
  // do not use hyperthreads
  int nthreads;
  ghost_machine_ncore(&nthreads,GHOST_NUMANODE_ANY);
  ghost_task_create(&task, nthreads, 0, void_lambda_caller, (void*) &context, 
        GHOST_TASK_DEFAULT, NULL, 0);
PHIST_DEB("enqueuing C++11-lambda as GHOST task and waiting for it\n");
  ghost_task_enqueue(task);
  ghost_task_wait(task);
  ghost_task_destroy(task);
}

// some helpful macros
/*
#define PHIST_GHOST_TASK_BEGIN phist_execute_lambda_as_ghost_task( [&]() {
#define PHIST_GHOST_TASK_END   } );
*/
#else /* PHIST_HAVE_CXX11_LAMBDAS */

#warning "C++11 not supported, not using GHOST tasking mechanism!"
#endif /* PHIST_HAVE_CXX11_LAMBDAS */

/* previously we started every kernel function as a task, but since this
   means that the threads need to be created every time, we now just start
   the main program as a task instead
*/   
#define PHIST_GHOST_TASK_BEGIN
#define PHIST_GHOST_TASK_END

#endif /* PHIST_GHOST_MACROS_HPP */
