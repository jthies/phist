#ifndef PHIST_KERNEL_FLAGS_H
#define PHIST_KERNEL_FLAGS_H

/* These flags allow the PHIST programmer to give the kernel lib
   "hints", but the kernel lib is not required to actually implement the feature.
   For instance, when creating a sparse matrix, the user can ask for a communication-
   reducing partitioning, but if the kernel lib does not support the feature, it may
   just ignore the flag.
   
   The integers 2^[0-15] are reserved for the kernel flags for now, but they may be aliased
   if there are no functions that accept two flags with the same integer value (e.g. a flag
   influencing cration of sparse matrices and a flag influencing the computation of a scalar
   product may both be defined as 32)
   
   In the documentation of the kernel functions, it is stated which flags are possible
   inputs to the particular function (combined by bitwise ors, '|').
   
   If no flags are documented for a function in in phist_kernels_decl.h, none are used by 
   any kernel lib.
*/
#define PHIST_IFLAG_DEFAULT 0

/* sparse matrix preprocessing: allow local matrix reorderings */
#define PHIST_SPARSEMAT_PERM_LOCAL 1
/* allow global symmetric permutations. If given, the partition sizes may  */
/* be adapted for load-balancing, and both global and local reordering may */
/* be used for e.g. reducing the bandwidth, number of messages and         */
/* communication volume                                                    */
#define PHIST_SPARSEMAT_PERM_GLOBAL 2
/* this flag is DEPRECATED, use SPARSEMAT_OPT_CARP instead */
/* in the current implementation of the builtin kernels,   */
/* DIST2_COLOR will cause a local permutation according to */
/* colors, which breaks the MPI communication but can be   */
/* used to assess the performance impact of coloring. A    */
/* warning will be issued.                                 */
#define PHIST_SPARSEMAT_DIST2_COLOR 4
#define PHIST_SPARSEMAT_OPT_SINGLESPMVM 8
#define PHIST_SPARSEMAT_OPT_BLOCKSPMVM 16
#define PHIST_SPARSEMAT_OPT_CARP 32
#define PHIST_SPARSEMAT_QUIET 64
/* if you pass pointers to maps to a function to create/read a matrix,  */
/* this flag makes the matrix responsible for deleting the map if it is */
/* deleted itself.                                                      */
/* If the creating function has no maps as arguments, they are created  */
/* and owned by the matrix automatically.                               */
#define PHIST_SPARSEMAT_OWN_MAPS 128

#define PHIST_SPARSEMAT_FLAGS_DESCRIPTION \
"     PHIST_IFLAG_DEFAULT 0 \n" \
"     PHIST_SPARSEMAT_PERM_LOCAL 1 \n" \
"     PHIST_SPARSEMAT_PERM_GLOBAL 2 \n" \
"     PHIST_SPARSEMAT_DIST2_COLOR 4 \n" \
"     PHIST_SPARSEMAT_OPT_SINGLESPMVM 8 \n" \
"     PHIST_SPARSEMAT_OPT_BLOCKSPMVM 16 \n" \
"     PHIST_SPARSEMAT_OPT_CARP 32 \n" \
"     PHIST_SPARSEMAT_QUIET 64 \n" \
"     PHIST_SPARSEMAT_OWN_MAPS 128 \n"

/* When this flag was passed to mvec_create, the memory for
   the multi-vector is allocated both on host and device for CUDA
   processes. This will allow using the functions mvec_extract_view,
   and mvec_from/to_device, which may return an error otherwise.
   
   Note that for sdMats we assume that the memory is allocated both on 
   host and device.
   */
#define PHIST_MVEC_REPLICATE_DEVICE_MEM 1

/* use more accurate reducitons or other floating point operations if available */
#define PHIST_ROBUST_REDUCTIONS 1

/* do not perform global MPI reduction on inner products */
#define PHIST_NO_GLOBAL_REDUCTION 2

/* sparseMat_times_mvec* flags, these are GHOST-specific up to now
   and should not be used in the code anywhere because they are subject
   to change. The purpose of these flags is benchmarking only.
*/
#define PHIST_SPMVM_ONLY_LOCAL 1024
#define PHIST_SPMVM_VECTOR 2048
#define PHIST_SPMVM_OVERLAP 4096
#define PHIST_SPMVM_TASK 8192

/* explicitly perform sdMat operation on the host CPU (the default). This flag affects kernels */
/* like sdMat_add_sdMat, sdMat_times_sdMat etc.                                                */
#define PHIST_SDMAT_RUN_ON_HOST 1024
/* execute sdMat function on device if this is a GPU process. This flag can e used to overrule  */
/* the default behavior of executing sdMat_add/times_sdMat and similar functions on the host.   */
/* The rationale behind not running sdMat operations on the GPU by default is that they usually */
/* involve few operations so the latency of launching a GPU kernel will be larger than doing    */
/* the computation on the host and uploading the result if it is needed afterwards on the GPU.  */
/*                                                                                              */
/* If PHIST_SDMAT_RUN_ON_HOST|PHIST_SDMAT_RUN_ON_DEVICE is speicified, the kernel is executed   */
/* on the host and the sdMat data is uploaded afterwards.                                       */
#define PHIST_SDMAT_RUN_ON_DEVICE 2048
/* shortcut for PHIST_SDMAT_RUN_ON_HOST|PHIST_SDMAT_RUN_ON_DEVICE */
#define PHIST_SDMAT_RUN_ON_HOST_AND_DEVICE (PHIST_SDMAT_RUN_ON_HOST|PHIST_SDMAT_RUN_ON_DEVICE)

#endif
