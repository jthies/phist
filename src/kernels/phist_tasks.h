#ifndef PHIST_TASKS_HPP
#define PHIST_TASKS_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"


/*! \def PHIST_TASK_DECLARE(taskName)
 *  Declare a tasks (e.g. a variable taskName of some type to be used later!)
 *  \taskName new identifier
 */

/*! \def PHIST_TASK_BEGIN(taskName)
 *  Marks the start of a task, must be followed by PHIST_TASK_END or PHIST_TASK_END_NOWAIT
 *  \param taskName previously declared identifier (with PHIST_TASK_DECLARE)
 */

/*! \def PHIST_TASK_END(task_ierr)
 * Marks the end of a task and waits for it to finish.
 * With GHOST this task allocates (and blocks) all available resources (e.g. cores)!
 * \warning currently only works in a void functions as it uses PHIST_CHK_IERR internally (so it may return with an error!)
 * \param tas_ierr (int*), for errors
 */

/*! \def PHIST_TASK_END_NOWAIT(task_ierr)
 * Marks the end of a task and runs it in the background
 * With GHOST this task has no resources (e.g. cores)
 * \warning currently only works in a void functions as it uses PHIST_CHK_IERR internally (so it may return with an error!)
 * \param tas_ierr (int*), for errors
 */

/*! \def PHIST_TASK_WAIT(taskName,task_ierr)
 * Waits for a background task to finish (e.g. one with PHIST_TASK_END_NOWAIT)
 * \warning currently only works in a void functions as it uses PHIST_CHK_IERR internally (so it may return with an error!)
 * \param taskName task identifier
 * \param tas_ierr (int*), for errors
 */

/*! \def PHIST_TASK_POST_STEP(task_ierr)
 * (Semaphore post) Marks some kind of progress in the *current* task (you can wait for it in *another* task with PHIST_TASK_WAIT_STEP)
 * \warning currently only works in a void functions as it uses PHIST_CHK_IERR internally (so it may return with an error!)
 * \param tas_ierr (int*), for errors
 */

/*! \def PHIST_TASK_WAIT_STEP(taskName,task_ierr)
 * (Semaphore wait) Wait for a PHIST_TASK_POST_STEP to be called in the specified task
 * \warning currently only works in a void functions as it uses PHIST_CHK_IERR internally (so it may return with an error!)
 * \param taskName task identifier
 * \param tas_ierr (int*), for errors
 */

/*! \def PHIST_MAIN_TASK_BEGIN
 * Starts the main task (if required). Put it directly after phist_kernels_init
 */

/*! \def PHIST_MAIN_TASK_END
 * Finishes the main task (if required). Put it directly before phist_kernels_finalize
 */


#if defined(__cplusplus) && defined(PHIST_HAVE_GHOST) && defined(PHIST_HAVE_CXX11_LAMBDAS)
#include <ghost.h>
#include <ghost/machine.h>


// glue between C++11 Lambda-functions and plain old C function pointers
namespace
{
template<typename LFunc>
class WrapLambdaForGhostTask
{
  public:
    // constructor
    WrapLambdaForGhostTask(ghost_task_t **task, const LFunc &context, int* iflag, bool async)
    {
      // check that the task is not already allocated/in use
      PHIST_CHK_IERR(*iflag = (*task != NULL) ? -1 : 0, *iflag);
      if( !async )
      {
        // don't create a thread at all if the current task is also not an async task
        ghost_task_t* curTask = NULL;
        ghost_task_cur(&curTask);
        if( curTask != NULL )
        {
          bool curTaskAsync = (curTask->nThreads == 0);
          if( !curTaskAsync )
          {
            // simply run the target code
            context();
            return;
          }
        }
      }


      // actually create the task and copy the context to the heap
      {
        int nThreads = async ? 0 : GHOST_TASK_FILL_ALL;
        ghost_task_flags_t flags = async ? GHOST_TASK_NOT_PIN : GHOST_TASK_DEFAULT;
        PHIST_CHK_GERR(ghost_task_create(task, nThreads, 0, void_lambda_caller, (void*) new LFunc(context), flags, NULL, 0), *iflag);
      }

      // If this task is not async, it needs resources.
      // If the parent has none, try to be adopted by the grand-parent, etc
      if( !async )
      {
        ghost_task_t *newParent = NULL;
        ghost_task_cur(&newParent);
        if( newParent )
        {
          while( newParent->nThreads == 0 )
          {
            if( newParent->parent == NULL )
              break;

            newParent = newParent->parent;
          }
        }
        (*task)->parent = newParent;
      }

      PHIST_DEB("enqueuing C++11-lambda as GHOST task and waiting for it\n");
      PHIST_CHK_GERR(ghost_task_enqueue(*task), *iflag);
      if( async )
        return;

      PHIST_CHK_GERR(ghost_task_wait(*task), *iflag);
      ghost_task_destroy(*task);
      *task = NULL;
    }

#if 0
    // C++11 move constructor
    WrapLambdaForGhostTask(WrapLambdaForGhostTask&&) = default;
#endif

  private:
    // hide default constructor etc
    WrapLambdaForGhostTask() = delete;
    WrapLambdaForGhostTask(const WrapLambdaForGhostTask&) = delete;
    WrapLambdaForGhostTask& operator=(const WrapLambdaForGhostTask&) = delete;

    // ghost requires a task function of the form "void* f(void*)"
    // a Lambda without capture can be converted automatically to a plain old C function pointer,
    // but the context goes out of scope when this function is left (so we copy it to the heap and destroy it here!)
    static void* void_lambda_caller(void* context)
    {
      (*(LFunc*)context)();               // run the target function
      delete (LFunc*)context;             // delete the copy of the lambda-function object
      return NULL;                        // we do not return anything
    }
};

// actually wrap a lambda for ghost and run it (creates an instance of WrapLambdaForGhostTask)
template<typename LFunc>
inline WrapLambdaForGhostTask<LFunc> phist_execute_lambda_as_ghost_task(ghost_task_t **task, const LFunc &context, int* iflag, bool async)
{
  return WrapLambdaForGhostTask<LFunc>(task,context,iflag,async);
}

// wait till a task is finished
inline void phist_wait_ghost_task(ghost_task_t** task, int* iflag)
{
  PHIST_CHK_IERR(*iflag = (*task != NULL) ? PHIST_SUCCESS : PHIST_WARNING, *iflag);

  // warn if we wait for a blocking task from within another blocking task -> deadlock
  ghost_task_t* curTask = NULL;
  ghost_task_cur(&curTask);
  if( curTask != NULL )
  {
    PHIST_CHK_IERR(*iflag = (*task)->nThreads > 0 && curTask->nThreads > 0 ? PHIST_WARNING : PHIST_SUCCESS, *iflag);
  }


  PHIST_CHK_GERR(ghost_task_wait(*task), *iflag);
  ghost_task_destroy(*task);
  *task = NULL;
}
}

// some helpful macros

#define PHIST_TASK_DECLARE(taskName) ghost_task_t* taskName = NULL;
#define PHIST_TASK_BEGIN(taskName) phist_execute_lambda_as_ghost_task(&taskName, ([&]() -> void {
#define PHIST_TASK_END(task_ierr)          }), task_ierr, false);PHIST_CHK_IERR((void)*(task_ierr),*(task_ierr));
#define PHIST_TASK_END_NOWAIT(task_ierr)   }), task_ierr, true );PHIST_CHK_IERR((void)*(task_ierr),*(task_ierr));
#define PHIST_TASK_WAIT(taskName,task_ierr) PHIST_CHK_IERR(phist_wait_ghost_task(&taskName,task_ierr),*(task_ierr));

#define PHIST_TASK_POST_STEP(task_ierr) {ghost_task_t* t = NULL; ghost_task_cur(&t); if( t != NULL ) sem_post(t->progressSem);}
#define PHIST_TASK_WAIT_STEP(taskName,task_ierr) {sem_wait(taskName->progressSem);}

// don't need a main task, every kernel function is put into a task!
#define PHIST_MAIN_TASK_BEGIN
#define PHIST_MAIN_TASK_END
//#define PHIST_MAIN_TASK_BEGIN {int task_ierr = 0; ghost_task_t* mainTask = NULL; phist_execute_lambda_as_ghost_task(&mainTask, [&]()->int {
//#define PHIST_MAIN_TASK_END   return 0;}, &task_ierr, false );PHIST_ICHK_IERR((void)task_ierr,task_ierr);}


#else /* __cplusplus && PHIST_HAVE_GHOST && PHIST_HAVE_CXX11_LAMBDAS */


#define PHIST_TASK_DECLARE(taskName)
#define PHIST_TASK_BEGIN(taskName)
#define PHIST_TASK_END(task_ierr)
#define PHIST_TASK_END_NOWAIT(task_ierr)
#define PHIST_TASK_WAIT(taskName,task_ierr)

#define PHIST_TASK_POST_STEP(task_ierr)
#define PHIST_TASK_WAIT_STEP(taskName,task_ierr)

#define PHIST_MAIN_TASK_BEGIN
#define PHIST_MAIN_TASK_END

#endif /* PHIST_HAVE_CXX11_LAMBDAS */

#endif /* PHIST_TASKS_HPP */
