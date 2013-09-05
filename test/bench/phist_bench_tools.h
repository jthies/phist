#ifdef PHIST_BENCH_TOOLS_H
#define PHIST_BENCH_TOOLS_H

#include "ghost.h"
#include "ghost_util.h"

#define PHIST_BENCH(stream,fcn_label,fcn) \
        { \
        double t0,t1; \
        const char* label="PHIST_BENCH"; \
        fprintf(stream,"%s FUNCTION %s\n",label,fcn_label);\
        fprintf(stream,"%s NUM_RANKS %d\n",label,ghost_getNumberOfRanks(MPI_COMM_WORLD));\
        fprintf(stream,"%s MAX_NUM_THREADS %d\n",label,ghost_thpool->nThreads);\
        fprintf(stream,"%s NUM_THREADS %d\n",label,ghost_ompGetNumThreads());\
        t0=ghost_wctime(); \
        fcn; \
        t1=ghost_wctime(); \
        fprintf(stream,"%s TIME %f\n",label,t1-t0);\
        }

#endif
