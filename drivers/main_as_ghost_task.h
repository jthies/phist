#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#ifdef PHIST_HAVE_LIKWID
#include "likwid.h"
#endif

#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#include "ghost/thpool.h"
#include "ghost/machine.h"
#include "ghost/pumap.h"
#include "ghost/locality.h"
#endif

#include "phist_macros.h"

#define PHIST_IPCHK_IERR(func,FLAG) {\
{func; if (FLAG!=PHIST_SUCCESS) { \
PHIST_OUT(PHIST_ERROR,"Error code %d (%s) returned from call %s\n(file %s, line %d)\n",\
(FLAG),(phist_retcode2str(FLAG)),(#func),(__FILE__),(__LINE__)); return &(FLAG);}}}

typedef struct mainArgs_t {
int argc;
char** argv;
} mainArgs_t;

static int* mainFunc(mainArgs_t* arg);

int main(int argc, char** argv)
{
  int ierr;
  mainArgs_t args;
  args.argc=argc;
  args.argv=argv;

#ifdef PHIST_HAVE_GHOST
  ghost_hwconfig_t hwconfig={.ncore=6,.nsmt=6};
  ghost_hwconfig_set(hwconfig);
#endif  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);
#ifdef PHIST_HAVE_GHOST
#ifndef PHIST_KERNEL_LIB_GHOST
    ghost_init(argc,argv);
#endif
    ghost_thpool_t *thpool;
    int nnuma = 0;
    int npu = 0;
    char *str;

    ghost_thpool_get(&thpool);
    ghost_machine_nnuma(&nnuma);
    ghost_machine_npu(&npu,GHOST_NUMANODE_ANY);

    ghost_pumap_string(&str);
    PHIST_OUT(PHIST_VERBOSE,"%s\n",str);
    free(str); str = NULL;

    PHIST_OUT(PHIST_VERBOSE,"The thread pool consists of %d threads\n",thpool->nThreads);

    ghost_task_t *mainTask;

    ghost_task_create(&mainTask, 10, 0, 
                (void*(*)(void*))&mainFunc, &args, GHOST_TASK_DEFAULT, NULL, 0);

    ghost_task_enqueue(mainTask);
    ghost_task_wait(mainTask);
    PHIST_ICHK_IERR(ierr=*((int*)mainTask->ret),ierr);
    ghost_task_destroy(mainTask);
#ifndef PHIST_KERNEL_LIB_GHOST
    ghost_finalize();
#endif
#else
// no ghost
ierr=*mainFunc(&args);
#endif    
  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
}

static int ierr;

