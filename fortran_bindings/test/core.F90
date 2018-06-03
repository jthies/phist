#include "phist_fort.h"
#include "phist_fort_test.h"

program phist_fort_core_test

use, intrinsic :: iso_c_binding
use mpi
use phist_types
use phist_tools_d
use phist_kernels
use phist_kernels_d
use phist_core_d

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

! matrix and vectors
TYPE(sparseMat_ptr) :: A
TYPE(mvec_ptr) :: x, y1, y2
! linear operator to wrap A
TYPE(DlinearOp), target :: A_op

integer(c_int), parameter :: ncols=8

! to compare the results we use the vector norms
real(c_double), dimension(ncols) :: norms1, norms2

real(c_double) :: alpha, beta, max_diff
integer :: rank
logical :: verbose

argc=0
argv=C_NULL_PTR

call MPI_Init(iflag)

! initialize test framework
call tests_init('core')

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

call phist_Dcreate_matrix(A,comm,C_CHAR_"graphene16x16"//C_NULL_CHAR,iflag)
ASSERT_EQ(0,iflag)


call phist_DsparseMat_get_domain_map(A,map,iflag)
ASSERT_EQ(0,iflag)
call phist_Dmvec_create(x,map,ncols,iflag)
ASSERT_EQ(0,iflag)
call phist_Dmvec_create(y1,map,ncols,iflag)
ASSERT_EQ(0,iflag)
call phist_Dmvec_create(y2,map,ncols,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_random(x,iflag)
ASSERT_EQ(0,iflag)
call phist_Dmvec_random(y1,iflag)
ASSERT_EQ(0,iflag)
! set y2=y1
call phist_Dmvec_add_mvec(1.D0,y1,0.D0,y2,iflag)
ASSERT_EQ(0,iflag)

call phist_DlinearOp_wrap_sparseMat(c_loc(A_op), A, iflag)
ASSERT_EQ(0,iflag)

! apply A and the wrapped A_op to x, and compare the norms of the output vectors
alpha=1.34834_c_double
beta=9.318583_c_double

call phist_DsparseMat_times_mvec(alpha,A,x,beta,y1,iflag)
EXPECT_EQ(0,iflag)

! currently there is no linearOp_apply function!!!
call phist_DlinearOp_apply(alpha,c_loc(A_op),x,beta,y2,iflag)
EXPECT_EQ(0,iflag)

! compute the norm of each column of y1 and y2 and compare them
call phist_Dmvec_norm2(y1,norms1(1),iflag);
EXPECT_EQ(0,iflag)

call phist_Dmvec_norm2(y2,norms2(1),iflag);
EXPECT_EQ(0,iflag)

max_diff=MAXVAL(ABS(norms1-norms2))
EXPECT_EQV(.true.,max_diff<1e-13)

call phist_DlinearOp_destroy(c_loc(A_op),iflag)
EXPECT_EQ(0,iflag)

call phist_DsparseMat_delete(A,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_delete(x,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_delete(y1,iflag)
ASSERT_EQ(0,iflag)

call phist_Dmvec_delete(y2,iflag)
ASSERT_EQ(0,iflag)

call tests_finalize(verbose)

call phist_kernels_finalize(iflag)

call MPI_Finalize(iflag)

end program phist_fort_core_test
