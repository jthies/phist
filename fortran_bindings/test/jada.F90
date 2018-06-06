#include "phist_fort.h"
#include "phist_fort_test.h"

program phist_fort_jada_test

use, intrinsic :: iso_c_binding
#ifdef PHIST_HAVE_MPI
use mpi
#endif
use phist_types
use phist_tools_d
use phist_kernels
use phist_kernels_d
use phist_core_d
use phist_jada
use phist_jada_d

use phist_testing

implicit none

! return flag
integer(c_int) :: iflag
! dummy args to phist_kernels_init
integer(c_int) :: argc
type(c_ptr), target :: argv

! communicator and map handles
TYPE(comm_ptr) :: comm
TYPE(const_map_ptr) :: map

! string for creating the matrix A
character(c_char), dimension(100) :: problem_string

! matrix and vectors
TYPE(sparseMat_ptr) :: A
TYPE(mvec_ptr) :: Q, v0
TYPE(sdMat_ptr) :: R


! linear operator to wrap A-sigma*I
! note: wrapped phist types are passed to phist functions
! by their C address (c_loc), which requires us to always
! declare them with the "target" qualifier. Alternatively
! we could just declare A_op of type(c_ptr) directly     
! since we don't really need access to its members here. 
TYPE(DlinearOp), target :: A_op


! options for the eigensolver
TYPE(jadaOpts), target :: opts
integer(c_int) :: nEigs, nEigsAlloc, nIter

complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: evals
real(c_double), allocatable, dimension(:)            :: resNorm

real(c_double) :: max_diff
integer :: rank
logical :: verbose

argc=0
argv=C_NULL_PTR
#ifdef PHIST_HAVE_MPI
call MPI_Init(iflag)
#endif
! initialize test framework
call tests_init('jada')

! it doesn't really matter what we pass in here
! *as long as we already initialized MPI!*
! Otherwise, the kernel lib may call the C variant
! of mpi_init and pass in these pointers
call phist_kernels_init(argc,c_loc(argv), iflag)
ASSERT_EQ(0,iflag)

! create a communicator object, note that this may be different from an MPI_Comm in Fortran
call phist_comm_create(comm,iflag)
ASSERT_EQ(0,iflag)

call phist_comm_get_rank(comm,rank,iflag)
ASSERT_EQ(0,iflag)

verbose = (rank==0)

call phist_Dcreate_matrix(A,comm,C_CHAR_"spinSZ16"//C_NULL_CHAR,iflag)
ASSERT_EQ(0,iflag)

call phist_jadaOpts_setDefaults(c_loc(opts))
ASSERT_EQ(0,iflag)

nEigs=10

opts%numEigs    = nEigs
opts%which      = phist_SR
opts%symmetry   = phist_HERMITIAN
opts%convTol        = 1.0e-6
opts%maxIters    = 30
opts%blockSize  = 2
opts%minBas     = 16
opts%maxBas    = 32
opts%innerSolvType      = phist_MINRES
opts%innerSolvMaxIters  = 10

! allocate memory for the solution
nEigsAlloc=nEigs+opts%blockSize-1
allocate(evals(nEigsAlloc), resNorm(nEigsAlloc))

call phist_DsparseMat_get_domain_map(A,map,iflag)
call phist_DsdMat_create(R,nEigsAlloc,nEigsAlloc,comm,iflag)
ASSERT_EQ(0,iflag)
call phist_Dmvec_create(Q,map,nEigsAlloc,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_create(v0,map,1,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_random(v0,iflag)
ASSERT_EQ(0,iflag)

opts%v0=v0

! we need a linearOp for A, but also - in contrast to the C interface - for B, even if it is the
! identity matrix. This is because we're passing the linearOp args by reference, so C_NULL_PTR  
! won't be accepted by the Fortran compiler (TODO)
call phist_DlinearOp_wrap_sparseMat(c_loc(A_op), A, iflag)
ASSERT_EQ(0,iflag)

! run the Jacobi-Davidson solver 'subspacejada'
call phist_Dsubspacejada(c_loc(A_op), C_NULL_PTR, opts,&
              Q, R, evals(1), resNorm(1), nEigs, nIter, iflag)              

! some assertions to test the interface worked
EXPECT_EQ(0,iflag)
EXPECT_EQV(.true.,nIter<opts%maxIters)
EXPECT_EQ(opts%numEigs,nEigs)

call phist_DlinearOp_destroy(c_loc(A_op),iflag)
ASSERT_EQ(0,iflag)

call phist_DsparseMat_delete(A,iflag)
ASSERT_EQ(0,iflag)

call phist_DsdMat_delete(R,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_delete(Q,iflag)
ASSERT_EQ(0,iflag)

call tests_finalize(verbose)

call phist_kernels_finalize(iflag)

#ifdef PHIST_HAVE_MPI
call MPI_Finalize(iflag)
#endif

end program phist_fort_jada_test
