#ifndef PHIST_KERNEL_FLAGS_H
#define PHIST_KERNEL_FLAGS_H

/* These flags allow the PHIST programmer to give the kernel lib
   "hints", but the kernel lib is not required to actually implement the feature.
   For instance, when creating a sparse matrix, the user can ask for a communication-
   reducing partitioning, but if the kernel lib does not support the feature, it may
   just ignore the flag.
   
   In the documentation of the kernel functions, it is stated which flags are possible
   inputs to the particular function (combined by bitwise ors, '|').
   
   If no flags are documented for a function in in phist_kernels_decl.h, none are used by 
   any kernel lib.
*/
#define PHIST_IFLAG_DEFAULT 0

/* sparse matrix preprocessing */
#define PHIST_SPARSEMAT_REPARTITION 1
#define PHIST_SPARSEMAT_DIST2_COLOR 2
#define PHIST_SPARSEMAT_OPT_SINGLESPMVM 4
#define PHIST_SPARSEMAT_OPT_BLOCKSPMVM 8

/* accuracy vs. speed */
#define PHIST_ROBUST_REDUCTIONS 1

#endif
