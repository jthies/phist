#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_macros.h"
#include <stdio.h>
#include <omp.h>

#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#include "ghost_taskq.h"
#else
#error "this driver makes no sense without ghost"
#endif

// this would be the 'block size' in JaDa
static const int ntasks=4;

typedef struct
  {
  int me;
  int num_threads;
  int core;
  int message;
  } arg_t;

static void *helloFunc(void* v_arg)
  {
  arg_t* arg = (arg_t*)v_arg;
#pragma omp parallel
  {
  omp_set_num_threads(ghost_thpool->nThreads);
  fprintf(stdout,"JOB %d enqueued by thread %d (of %d) on core %d\n run by thread %d (of %d) on core %d\n",
        arg->message,
        arg->me,arg->num_threads,arg->core,
        omp_get_thread_num(),omp_get_num_threads(),ghost_getCore());  
  }
  return NULL;
  }

static void Hello_World_via_the_Queue(int message)
  {
  arg_t arg;
  ghost_task_t *task;
  arg.me = omp_get_thread_num();
  arg.core=ghost_getCore();
  arg.message=message;
  arg.num_threads = omp_get_num_threads();
  task = ghost_task_init(GHOST_TASK_FILL_ALL, 0, &helloFunc, (void*)&arg, 
  GHOST_TASK_DEFAULT);

  ghost_task_add(task);

  ghost_task_wait(task);
  ghost_task_destroy(task);
  }

  
int main(int argc, char** argv)
  {
  int rank, num_proc;
  int ierr;
  int i;

  comm_ptr_t comm_world;
  
  PHIST_CHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);
  // avoid duplicate init call
#ifndef PHIST_KERNEL_LIB_GHOST
  // initialize ghost queue
  ghost_init(argc,argv);
#endif
  PHIST_CHK_IERR(phist_comm_create(&comm_world,&ierr),ierr);
  
  fprintf(stdout,"num threads seen by ghost: %d\n",ghost_thpool->nThreads);
  fprintf(stdout,"max num threads seen by OpenMP: %d\n",omp_get_max_threads());

omp_set_num_threads(ntasks);

// every task thread starts one instance of a function,
// and executes it right away:
#pragma omp parallel private(i)
  {
  i=omp_get_thread_num();
  fprintf(stdout,"control thread %d of %d starting a task\n",i,omp_get_num_threads());
#pragma omp task
  Hello_World_via_the_Queue(i);
  }
  PHIST_CHK_IERR(phist_kernels_finalize(&ierr),ierr);
  // avoid duplicate finish call
#ifndef PHIST_KERNEL_LIB_GHOST
  // finalize ghost queue
  ghost_finish();
#endif  
  }
