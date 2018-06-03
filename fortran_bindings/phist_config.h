#ifndef PHIST_CONFIG_H
#define PHIST_CONFIG_H

/* access the version of phist as a string: "major.minor.patch" */
#define PHIST_VERSION_STRING "1.6.1"

/* components of the version number */
#define PHIST_VERSION_MAJOR 1
#define PHIST_VERSION_MINOR 6
#define PHIST_VERSION_PATCH  1

/* access the version of phist: major|minor|patch */
#define PHIST_VERSION_INT 10601

/* git revision at which cmake was invoked */
#define PHIST_GIT_REVISION "b9f132e" 

#define PHIST_INSTALL_INFO \
"CMAKE_VERSION: 3.10.1\n" \
"CXX_COMPILER: /home/jonas/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/gcc-7.3.0-5njqlxwe45dwe6xqz6rsfcalihnumjh5/bin/gcc\n" \
"CXX_COMPILER: /home/jonas/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/gcc-7.3.0-5njqlxwe45dwe6xqz6rsfcalihnumjh5/bin/g++\n" \
"Fortran_COMPILER: /home/jonas/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/gcc-7.3.0-5njqlxwe45dwe6xqz6rsfcalihnumjh5/bin/gfortran\n" \
"CXX_FLAGS:      -fopenmp -fopenmp-simd\n" \
"C_FLAGS:        -fopenmp -fopenmp-simd -msse -mavx\n" \
"Fortran_FLAGS: -cpp -ffree-line-length-none -fopenmp -fopenmp-simd\n"

/* general verbosity of PHIST:
        0: only print errors
        1: errors and warnings
        2: some info about the progress of the algorthms etc.
        3: verbose info
        4: debugging output
        5: function tracing
*/
#define PHIST_OUTLEV 2

/* enables some potentially expensive tests, for instance whether the index spaces of
   input and output vectors to an operation are correct. By default, this is ON if
   CMAKE_BUILD_TYPE=Debug.
 */
/* #undef PHIST_TESTING */

/* some compiler/system flags */
#define PHIST_HAVE_CXX11_LAMBDAS
#define PHIST_HAVE_CXX11_THREADLOCAL
#define PHIST_HAVE_CXX11_MOVEDEFAULT
#define PHIST_HAVE_OPENMP_SIMD
#define PHIST_HAVE_MPI
#define PHIST_HAVE_OPENMP

/* cross-compile for Intel MIC */
/* #undef PHIST_BUILD_MIC */

/* do you want to compile the single-precision
version of the code in addition to the standard 
(double)? 
*/
/* #undef PHIST_HAVE_SP */

/* by default PHIST tries to use 64-bit indices with Epetra/Tpetra.
   This option can be used to force 32-bit (int) global indices, this
   may be useful for interfaceing e.g. with partitioning libraries like
   Zoltan/Isorropia. This parameter has no effect on the builtin and ghost
   kernel interface.
 */
/* #undef PHIST_FORCE_32BIT_GIDX */

/* do you want to compile the complex valued
version of the code (only available for some
kernel libs) ?
*/
/* #undef PHIST_HAVE_CMPLX */

/* one can enable/disable the use of certain SIMD intrinsics (manually) using ccmake,
 * this mostly affects the builtin kernels, if you e.g. use GHOST as kernel lib, it
 * will use whatever it was configured to use when installed.
 */
#define PHIST_HAVE_SSE
#define PHIST_HAVE_AVX
/* #undef PHIST_HAVE_AVX2 */
/* #undef PHIST_HAVE_AVX512 */

/* do we have support for higher precision kernel routines (Kahan-style reductions etc.)?
*/
/* #undef PHIST_HIGH_PRECISION_KERNELS */
/* #undef PHIST_HIGH_PRECISION_KERNELS_FORCE */

/* which library provides the basic matrix/vector operations? */
/* #undef PHIST_KERNEL_LIB_EPETRA */
/* #undef PHIST_KERNEL_LIB_TPETRA */
/* #undef PHIST_KERNEL_LIB_GHOST */
#define PHIST_KERNEL_LIB_BUILTIN
/* #undef PHIST_KERNEL_LIB_MAGMA */
/* #undef PHIST_KERNEL_LIB_PETSC */
/* #undef PHIST_KERNEL_LIB_EIGEN */

#define PHIST_MVECS_ROW_MAJOR

/* use random number generator prescribed by the kernel lib or our own? 
   The builtin variant may be much faster but the PHIST RNG gives reproducible
   results across runs with different number of processes/threads 
   */
#define PHIST_BUILTIN_RNG

/* ESSEX libraries: optional */
/* #undef PHIST_HAVE_GHOST */
/* #undef PHIST_HAVE_ESSEX_PHYSICS */

/* third party libraries: optional */

/* Trilinos */
/* #undef PHIST_HAVE_TEUCHOS */
/* #undef PHIST_HAVE_KOKKOS */
/* #undef PHIST_HAVE_ANASAZI */
/* #undef PHIST_HAVE_BELOS */
/* #undef PHIST_HAVE_IFPACK */
/* #undef PHIST_HAVE_ML */
/* #undef PHIST_HAVE_IFPACK2 */
/* #undef PHIST_HAVE_MUELU */
/* #undef PHIST_HAVE_AMESOS2 */

/* #undef PHIST_HAVE_ISORROPIA */
/* #undef PHIST_HAVE_ZOLTAN */

/* other TPLs */
/* #undef PHIST_HAVE_MKL */
/* #undef PHIST_HAVE_PARMETIS */
/* #undef PHIST_HAVE_KAHIP */
/* #undef PHIST_HAVE_COLPACK */
/* #undef PHIST_HAVE_LIKWID */
/* do we have access to MPACK (multi-precision BLAS/LAPACK and GMP (GNU Multiple Precision library)
 * This is currently not used/implemented.
 */
/* #undef PHIST_HAVE_MPACK_GMP */
/* MPACK with QD library, used in robust SVQB with builtin kernels 
 */
/* #undef PHIST_HAVE_MPACK_QD */

/* TPL-specific settings */
/* #undef PHIST_ANASAZI_THEIR_ORTHO_MANAGER */
/* #undef PHIST_BELOS_THEIR_ORTHO_MANAGER */

/* timing/profiling/benchmarking behavior */
/* #undef LIKWID_PERFMON */
#define PHIST_TIMEMONITOR
/* #undef PHIST_TIMEMONITOR_PERLINE */
/* #undef PHIST_USE_TEUCHOS_TIMEMONITOR */
/* #undef PHIST_TIMINGS_FULL_TRACE */
/* #undef PHIST_PERFCHECK */
/* #undef PHIST_PERFCHECK_REALISTIC */
/* #undef PHIST_PERFCHECK_SEPARATE_OUTPUT */
#define PHIST_BENCH_LARGE_N 134217728
#define PHIST_PEAK_FLOPS_CPU -1
#define PHIST_PEAK_FLOPS_GPU -1

/* specific settings to tune the GHOST performance */
#ifdef PHIST_KERNEL_LIB_GHOST
/* with GHOST we execute everything in tasks, this can make debugging a bit tedious
   because GDB will jump out of the function context and into the tasking stuff. If
   this flag is defined, all kernels are executed in order.
 */
/* #undef PHIST_USE_GHOST_TASKS */
/* #undef PHIST_SELL_C */
/* #undef PHIST_SELL_SIGMA */
/* default mode flag passed to ghost_spmv, see ghost/spmv.h for details.
   as of GHOST 1.1:
   0: vector mode
   8: overlap mode (using non-blocking MPI)
  16: task mode (using ghost tasks for the overlapping)
  
  The default value can be overridden by flags to phist_XsparseMat_times_mvec*, see
  src/kernels/phist_kernel_flags.h
*/
#define PHIST_DEFAULT_SPMV_MODE 
#endif

#define PHIST_TRY_TO_PIN_THREADS

# ifdef PHIST_HAVE_GHOST
# include "ghost/config.h"
#  if defined(GHOST_HAVE_MPI)&&(!defined(PHIST_HAVE_MPI))
#  error "ghost was compiled with MPI, conflict with phist installation without MPI."
#  endif
#  if (!defined(GHOST_HAVE_MPI))&&defined(PHIST_HAVE_MPI)
#  error "ghost was compiled without MPI, conflict with phist installation with MPI."
#  endif
# endif

#endif
