#ifndef PHIST_GHOST_MACROS_HPP
#define PHIST_GHOST_MACROS_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <ghost.h>


#ifdef PHIST_HAVE_CXX11_LAMBDAS

// glue between C++11 Lambda-functions and plain old C function pointers
template<typename LFunc>
void phist_execute_lambda_as_ghost_task(LFunc context)
{
  // a Lambda without capture can be converted automatically to a plain old C function pointer
  auto void_lambda_caller = [] (void* context) -> void* {(*(LFunc*)context)(); return NULL;};

  ghost_task_t *task = NULL;
  ghost_task_create(&task, GHOST_TASK_FILL_ALL, 0, void_lambda_caller, (void*) &context, GHOST_TASK_DEFAULT, NULL, 0);
PHIST_DEB("enqueuing C++11-lambda as GHOST task and waiting for it\n");
  ghost_task_enqueue(task);
  ghost_task_wait(task);
  ghost_task_destroy(task);
}

// some helpful macros
#define PHIST_GHOST_TASK_BEGIN phist_execute_lambda_as_ghost_task( [&]() {
#define PHIST_GHOST_TASK_END   } );

#else /* PHIST_HAVE_CXX11_LAMBDAS */

#warning "C++11 not supported, not using GHOST tasking mechanism!"
#define PHIST_GHOST_TASK_BEGIN
#define PHIST_GHOST_TASK_END

#endif /* PHIST_HAVE_CXX11_LAMBDAS */


#endif /* PHIST_GHOST_MACROS_HPP */
