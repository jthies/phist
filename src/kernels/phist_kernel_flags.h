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
/* this flag is DEPRECATED, use SPARSEMAT_OPT_CARP instead */
/* in the current implementation of the builtin kernels,   */
/* DIST2_COLOR will cause a local permutation according to */
/* colors, which breaks the MPI communication but can be   */
/* used to assess the performance impact of coloring. A    */
/* warning will be issued.                                 */
#define PHIST_SPARSEMAT_DIST2_COLOR 2
#define PHIST_SPARSEMAT_OPT_SINGLESPMVM 4
#define PHIST_SPARSEMAT_OPT_BLOCKSPMVM 8
#define PHIST_SPARSEMAT_OPT_CARP 16
#define PHIST_SPARSEMAT_QUIET 32

#define PHIST_SPARSEMAT_FLAGS_DESCRIPTION \
"     PHIST_IFLAG_DEFAULT 0 \n" \
"     PHIST_SPARSEMAT_REPARTITION 1 \n" \
"     PHIST_SPARSEMAT_DIST2_COLOR 2 \n" \
"     PHIST_SPARSEMAT_OPT_SINGLESPMVM 4 \n" \
"     PHIST_SPARSEMAT_OPT_BLOCKSPMVM 8 \n" \
"     PHIST_SPARSEMAT_OPT_CARP 16 \n" \
"     PHIST_SPARSEMAT_QUIET 32 \n"

/* When this flag was passed to mvec_create, the memory for
   the multi-vector is allocated both on host and device for CUDA
   processes. This will allow using the functions mvec_extract_view,
   and mvec_from/to_device, which may return an error otherwise.
   
   Note that for sdMats we assume that the memory is allocated both on 
   host and device.
   */
#define PHIST_MVEC_REPLICATE_DEVICE_MEM 1

/* sparseMat_times_mvec* flags, these are GHOST-sepcific up to now
   and should not be used in the code anywhere because they are subject to
   to change. The purpose of these flags is benchmarking only.
   */
#define PHIST_SPMVM_ONLY_LOCAL 2
#define PHIST_SPMVM_VECTOR 4
#define PHIST_SPMVM_OVERLAP 8
#define PHIST_SPMVM_TASK 16

/* accuracy vs. speed */
#define PHIST_ROBUST_REDUCTIONS 1

#endif
