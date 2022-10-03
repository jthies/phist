#include "phist_fort.h"
#include "phist_fort_test.h"

program phist_fort_kernels_test

use, intrinsic :: iso_c_binding
#ifdef PHIST_HAVE_MPI
use mpi
#endif
use phist_types
use phist_kernels
use phist_kernels_s
use phist_kernels_d
use phist_kernels_c
use phist_kernels_z

use phist_testing

implicit none

! return flag
integer(c_int) :: iflag
! dummy args to phist_kernels_init
integer(c_int) :: argc
type(c_ptr), target, dimension(1) :: argv

! indicate which types are available in the kernel library
logical :: haveS, haveC, haveD, haveZ

! communicator handle
TYPE(comm_ptr) :: comm

! for testing sdMat functionality
type(sdMat_ptr) :: M1, M2
type(sdMat_ptr) :: M_view
integer(lidx) :: nrowsM, ncolsM
integer(lidx) :: nrowsM_view, ncolsM_view
integer(lidx) :: ldM, ldM_view

integer(lidx) :: i0, i1, j0, j1
integer, dimension(2) :: shape

type(c_ptr) :: M_dat_C, M_view_dat_C
real(c_double), pointer, dimension(:,:) :: M_dat_F, M_View_dat_F

real(c_double) :: max_diff
integer :: rank
logical :: verbose

argc=0
argv(1)=C_NULL_PTR

#ifdef PHIST_HAVE_MPI
call MPI_Init(iflag)
#endif
! initialize test framework
call tests_init('kernels')

! it doesn't really matter what we pass in here
! *as long as we already initialized MPI!*
! Otherwise, the kernel lib may call the C variant
! of mpi_init and pass in these pointers
call phist_kernels_init(argc,argv, iflag)
ASSERT_EQ(0,iflag)

! create a communicator object, note that this may be different from an MPI_Comm in Fortran
call phist_comm_create(comm,iflag)
ASSERT_EQ(0,iflag)

call phist_comm_get_rank(comm,rank,iflag)
ASSERT_EQ(0,iflag)

verbose = (rank==0)

! check which data types are supported by the phist installation
#ifdef PHIST_HAVE_SP
call phist_Stype_avail(iflag);
haveS=(iflag==0);
if (haveS) then
  write(*,*) 'have S (float)'
end if
EXPECT_EQV(.true.,haveS)
#endif

call phist_Dtype_avail(iflag);
haveD=(iflag==0);
if (haveD) then
  write(*,*) 'have D (double)'
end if
EXPECT_EQV(.true.,haveD)

#if defined(PHIST_HAVE_SP)&&defined(PHIST_HAVE_CMPLX)
call phist_Ctype_avail(iflag);
haveC=(iflag==0);
if (haveC) then
  write(*,*) 'have C (complex float)'
end if
EXPECT_EQV(.true.,haveC)
#endif
#ifdef PHIST_HAVE_CMPLX
call phist_Ztype_avail(iflag);
haveZ=(iflag==0);
if (haveZ) then
  write(*,*) 'have Z (complex double)'
end if
EXPECT_EQV(.true.,haveZ)
#endif

! sdMat tests
nrowsM=5
ncolsM=6
call phist_DsdMat_create(M1,nrowsM,ncolsM,comm,iflag)
ASSERT_EQ(0,iflag)

call phist_DsdMat_create(M2,nrowsM,ncolsM,comm,iflag)
ASSERT_EQ(0,iflag)

M_view=C_NULL_PTR
! NOTE: C-like 0-based indexing, so viewing columns[i0:i1,j0:j1]
!       skips the first and last row and column if:
i0=1
i1=nrowsM-2
j0=1
j1=ncolsM-2
nrowsM_view=i1-i0+1
ncolsM_view=i1-i0+1
call phist_DsdMat_view_block(M1, M_view, i0, i1, j0, j1, iflag)
ASSERT_EQ(0,iflag)

! get pointers to actual entries, this requires an extra step to
! tell Fortran a out the shape of the data
call phist_DsdMat_extract_view(M1,M_dat_C,ldM,iflag)
ASSERT_EQ(0,iflag)
shape(1)=ldM
shape(2)=ncolsM
call c_f_pointer(M_dat_C,M_dat_F,shape)

call phist_DsdMat_extract_view(M_view,M_view_dat_C,ldM_view,iflag)
ASSERT_EQ(0,iflag)
ASSERT_EQ(ldM,ldM_view)

shape(1)=ldM_view
shape(2)=ncolsM_view
call c_f_pointer(M_view_dat_C,M_view_dat_F,shape)

! set the complete matrix to 1
call phist_DsdMat_put_value(M1,1.0_8,iflag)

! set the viewed part to 2
call phist_DsdMat_put_value(M_view,2.0_8,iflag)

if (verbose) then
        write(*,*) 'sdMat data in Fortran:'
        write(*,'(6F6.1)') transpose(M_dat_F(1:nrowsM,1:ncolsM))
end if

max_diff=MAXVAL(ABS(M_dat_F(i0+1:i1+1,j0+1:j1+1)-2.0))
EXPECT_EQ(0.0,max_diff)
max_diff=MAXVAL(ABS(M_view_dat_F(1:nrowsM_view,1:ncolsM_view)-2.0))
EXPECT_EQ(0.0,max_diff)
max_diff=MAXVAL(ABS(M_dat_F(1:nrowsM,1)-1.0))
EXPECT_EQ(0.0,max_diff)
max_diff=MAXVAL(ABS(M_dat_F(1:nrowsM,ncolsM)-1.0))
EXPECT_EQ(0.0,max_diff)
max_diff=MAXVAL(ABS(M_dat_F(1,1:ncolsM)-1.0))
EXPECT_EQ(0.0,max_diff)
max_diff=MAXVAL(ABS(M_dat_F(nrowsM,1:ncolsM)-1.0))
EXPECT_EQ(0.0,max_diff)

call phist_DsdMat_delete(M_view,iflag)
ASSERT_EQ(0,iflag)

call phist_DsdMat_delete(M1,iflag)
ASSERT_EQ(0,iflag)

call phist_DsdMat_delete(M2,iflag)
ASSERT_EQ(0,iflag)

call tests_finalize(verbose)

call phist_kernels_finalize(iflag)

#ifdef PHIST_HAVE_MPI
call MPI_Finalize(iflag)
#endif
end program phist_fort_kernels_test
