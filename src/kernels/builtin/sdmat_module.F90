/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

!> \file sdmat_module.f90
!! Defines sdmat_module, the phist builtin implementation of phist_DsdMat_*
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

#include "phist_config_fortran.h"

!> Implementation of phist_DsdMat_* of the phist builtin kernels
!!
!! Everything should be straight-forward - not tuned specifically for
!! performance as it is assumed that sdMats are much smaller mvecs
module sdmat_module
  use mpi_f08, only: MPI_Comm
  implicit none
  private

  public :: SDMat_t
  !public :: phist_DsdMat_create
  !public :: phist_DsdMat_delete
  !public :: phist_DsdMat_extract_view
  !public :: phist_DsdMat_get_nrows
  !public :: phist_DsdMat_get_ncols
  !public :: phist_DsdMat_view_block
  !public :: phist_DsdMat_get_block
  !public :: phist_DsdMat_set_block
  !public :: phist_DsdMat_put_value
  !public :: phist_DsdMat_print
  !public :: phist_DsdMat_random
  !public :: phist_DsdMat_identity
  !public :: phist_DsdMat_add_sdMat
  public :: sdmat_add_sdmat
  !public :: phist_DsdMat_times_sdMat
  !public :: phist_DsdMatT_times_sdMat
  !public :: phist_DsdMat_times_sdMatT
  public :: sdmat_times_sdmat


  !==================================================================================
  !> sdmat with col-wise layout
  !! also stores the error for more precise operations
  type SDMat_t
    !--------------------------------------------------------------------------------
    integer     :: imin, imax
    integer     :: jmin, jmax
    type(MPI_Comm) :: comm
    real(kind=8), contiguous, pointer :: val(:,:) => null()
#ifdef PHIST_HIGH_PRECISION_KERNELS
    real(kind=8), contiguous, pointer :: err(:,:) => null()
#endif
    logical     :: is_view
    !--------------------------------------------------------------------------------
  end type SDMat_t


  !==================================================================================
  ! required interfaces of C functions
  interface
    !void phist_Drandom_number(int n, double*x);
    subroutine phist_Drandom_number(n,x) bind(C,name='phist_Drandom_number')
      use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), intent(out) :: x
    end subroutine
#ifdef PHIST_HIGH_PRECISION_KERNELS
    !void daxpby_prec(int n, double alpha, const double *restrict a, const double *restrict aC,
    !                        double beta,        double *restrict b,       double *restrict bC)
    subroutine daxpby_prec(n,alpha,a,aC,beta,b,bC) bind(C)
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value :: n
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: a(*), aC(*)
      real(kind=C_DOUBLE), intent(inout) :: b(*), bC(*)
    end subroutine
    !void dgemm_prec(int m, int n, int k, double alpha, const double *restrict a, const double *restrict aC,
    !                                                   const double *restrict b, const double *restrict bC,
    !                                     double beta,        double *restrict c,       double *restrict cC)
    subroutine dgemm_prec(m,n,k, alpha,a,aC,b,bC, beta,c,cC) bind(C)
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value :: m,n,k
      real(kind=C_DOUBLE), value :: alpha, beta
      real(kind=C_DOUBLE), intent(in) :: a(*),aC(*),b(*),bC(*)
      real(kind=C_DOUBLE), intent(inout) :: c(*),cC(*)
    end subroutine
#endif
  end interface
contains

  !==================================================================================
  ! axpby operation for SDMat_t
  subroutine sdmat_add_sdmat(alpha,A,beta,B)
    !--------------------------------------------------------------------------------
    real(kind=8),   intent(in)    :: alpha
    type(SDMat_t),  intent(in)    :: A
    real(kind=8),   intent(in)    :: beta
    type(SDMat_t),  intent(inout) :: B
    !--------------------------------------------------------------------------------
    real(kind=8), allocatable :: a_(:,:), b_(:,:), aC_(:,:), bC_(:,:)
    integer :: nr, nc
    !--------------------------------------------------------------------------------

    nr = B%imax-B%imin+1
    nc = B%jmax-B%jmin+1
    allocate(a_(nr,nc),b_(nr,nc),aC_(nr,nc),bC_(nr,nc))
    a_  = A%val(A%imin:A%imax,A%jmin:A%jmax)
    b_  = B%val(B%imin:B%imax,B%jmin:B%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    aC_ = A%err(A%imin:A%imax,A%jmin:A%jmax)
    bC_ = B%err(B%imin:B%imax,B%jmin:B%jmax)
#endif

#ifdef PHIST_HIGH_PRECISION_KERNELS
    call daxpby_prec(nr*nc,alpha,a_,aC_,beta,b_,bC_)
#else
    b_ = alpha*a_+beta*b_
#endif

    B%val(B%imin:B%imax,B%jmin:B%jmax) = b_
#ifdef PHIST_HIGH_PRECISION_KERNELS
    B%err(B%imin:B%imax,B%jmin:B%jmax) = bC_
#endif

  end subroutine sdmat_add_sdmat


  !==================================================================================
  ! axpby operation for SDMat_t
  subroutine sdmatT_add_sdmat(alpha,A,beta,B)
    !--------------------------------------------------------------------------------
    real(kind=8),   intent(in)    :: alpha
    type(SDMat_t),  intent(in)    :: A
    real(kind=8),   intent(in)    :: beta
    type(SDMat_t),  intent(inout) :: B
    !--------------------------------------------------------------------------------
    real(kind=8), allocatable :: a_(:,:), b_(:,:), aC_(:,:), bC_(:,:)
    integer :: nr, nc
    !--------------------------------------------------------------------------------

    nr = B%imax-B%imin+1
    nc = B%jmax-B%jmin+1
    allocate(a_(nr,nc),b_(nr,nc),aC_(nr,nc),bC_(nr,nc))
    a_  = transpose(A%val(A%imin:A%imax,A%jmin:A%jmax))
    b_  = B%val(B%imin:B%imax,B%jmin:B%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    aC_ = transpose(A%err(A%imin:A%imax,A%jmin:A%jmax))
    bC_ = B%err(B%imin:B%imax,B%jmin:B%jmax)
#endif

#ifdef PHIST_HIGH_PRECISION_KERNELS
    call daxpby_prec(nr*nc,alpha,a_,aC_,beta,b_,bC_)
#else
    b_ = alpha*a_+beta*b_
#endif

    B%val(B%imin:B%imax,B%jmin:B%jmax) = b_
#ifdef PHIST_HIGH_PRECISION_KERNELS
    B%err(B%imin:B%imax,B%jmin:B%jmax) = bC_
#endif

  end subroutine sdmatT_add_sdmat


  !==================================================================================
  ! gemm operation for SDMat_t
  subroutine sdmat_times_sdmat(transA,transB,alpha,A,B,beta,C,ierr)
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    character(len=1),intent(in)   :: transA, transB
    real(kind=8),   intent(in)    :: alpha
    type(SDMat_t),  intent(in)    :: A
    type(SDMat_t),  intent(in)    :: B
    real(kind=8),   intent(in)    :: beta
    type(SDMat_t),  intent(inout) :: C
    integer(kind=C_INT), intent(out)  :: ierr
    !--------------------------------------------------------------------------------
    integer :: m, n, k
    real(kind=8), allocatable :: a_(:,:), b_(:,:), aC_(:,:), bC_(:,:), c_(:,:), cC_(:,:)
    !--------------------------------------------------------------------------------

    ! determine data layout
    m = C%imax-C%imin+1
    n = C%jmax-C%jmin+1
    if( transA .eq. 'T' .or. transA .eq. 't') then
      k = A%imax-A%imin+1
    else
      k = A%jmax-A%jmin+1
    end if

    allocate(a_(m,k),b_(k,n),c_(m,n),aC_(m,k),bC_(k,n),cC_(m,n))
    if( transA .eq. 'T' .or. transA .eq. 't') then
      a_  = transpose(A%val(A%imin:A%imax,A%jmin:A%jmax))
#ifdef PHIST_HIGH_PRECISION_KERNELS
      aC_ = transpose(A%err(A%imin:A%imax,A%jmin:A%jmax))
#endif
    else
      a_  = A%val(A%imin:A%imax,A%jmin:A%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
      aC_ = A%err(A%imin:A%imax,A%jmin:A%jmax)
#endif
    end if
    if( transB .eq. 'T' .or. transB .eq. 't') then
      b_  = transpose(B%val(B%imin:B%imax,B%jmin:B%jmax))
#ifdef PHIST_HIGH_PRECISION_KERNELS
      bC_ = transpose(B%err(B%imin:B%imax,B%jmin:B%jmax))
#endif
    else
      b_  = B%val(B%imin:B%imax,B%jmin:B%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
      bC_ = B%err(B%imin:B%imax,B%jmin:B%jmax)
#endif
    end if
    c_  = C%val(C%imin:C%imax,C%jmin:C%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    cC_ = C%err(C%imin:C%imax,C%jmin:C%jmax)
#endif

#ifdef PHIST_HIGH_PRECISION_KERNELS
    call dgemm_prec(m,n,k, alpha,a_,aC_,b_,bC_, beta,c_,cC_)
#else
    c_ = alpha*matmul(a_,b_) + beta*c_
#endif

    C%val(C%imin:C%imax,C%jmin:C%jmax) = c_
#ifdef PHIST_HIGH_PRECISION_KERNELS
    C%err(C%imin:C%imax,C%jmin:C%jmax) = cC_
#endif

    ierr = 0

  end subroutine sdmat_times_sdmat

  !==================================================================================
  ! wrapper routines 

  subroutine phist_DsdMat_create(sdmat_ptr, nrows, ncols, comm_ptr, ierr) bind(C,name='phist_DsdMat_create_f')
    use, intrinsic :: iso_c_binding
    use mpi_f08
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: sdmat_ptr
    integer(C_INT),     value       :: nrows, ncols
    type(C_PTR),        value       :: comm_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    type(MPI_Comm), pointer :: comm
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

    allocate(sdmat)
    sdmat%is_view = .false.
    sdmat%imin = 1
    sdmat%imax = nrows
    sdmat%jmin = 1
    sdmat%jmax = ncols
    if( c_associated(comm_ptr) ) then
      call c_f_pointer(comm_ptr, comm)
      sdmat%comm = comm
    else
      sdmat%comm = MPI_COMM_SELF
    end if
#ifdef F_DEBUG
    write(*,*) 'creating new sdmat with dimensions:', nrows, ncols, 'address', transfer(c_loc(sdmat),dummy)
    flush(6)
#endif
    allocate(sdmat%val(nrows,ncols))
    sdmat%val = 0._8
#ifdef PHIST_HIGH_PRECISION_KERNELS
    allocate(sdmat%err(nrows,ncols))
    sdmat%err = 0._8
#endif
    sdmat_ptr = c_loc(sdmat)
    ierr = 0

  end subroutine phist_DsdMat_create

  subroutine phist_DsdMat_create_view(sdmat_ptr, comm_ptr,c_val,lda,nrows,ncols,ierr) &
  bind(C,name='phist_DsdMat_create_view_f')
    use, intrinsic :: iso_c_binding
    use mpi_f08
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: sdmat_ptr
    type(C_PTR),        value       :: comm_ptr
    integer(C_INT64_T), value       :: lda
    integer(C_INT),     value       :: nrows, ncols
    REAL(c_double), target          :: c_val(lda,ncols)
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    type(MPI_Comm), pointer :: comm
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

    allocate(sdmat)
    sdmat%is_view = .true.
    
    sdmat%imin = 1
    sdmat%imax = nrows
    sdmat%jmin = 1
    sdmat%jmax = ncols
    if( c_associated(comm_ptr) ) then
      call c_f_pointer(comm_ptr, comm)
      sdmat%comm = comm
    else
      sdmat%comm = MPI_COMM_NULL
    end if
#ifdef F_DEBUG
    write(*,*) 'creating sdmat view with dimensions:', nrows, ncols, 'address', transfer(c_loc(sdmat),dummy)
    flush(6)
#endif
    sdmat%val=>c_val
#ifdef PHIST_HIGH_PRECISION_KERNELS
    allocate(sdmat%err(nrows,ncols))
#endif
    sdmat_ptr = c_loc(sdmat)
    ierr = 0

  end subroutine phist_DsdMat_create_view


  subroutine phist_DsdMat_delete(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_delete_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: sdmat_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

#ifdef F_DEBUG
    write(*,*) 'deleting sdmat at address', transfer(sdmat_ptr,dummy)
    flush(6)
#endif
    if( c_associated(sdmat_ptr) ) then
      call c_f_pointer(sdmat_ptr, sdmat)
      if( .not. sdmat%is_view) then
        deallocate(sdmat%val)
#ifdef PHIST_HIGH_PRECISION_KERNELS
        deallocate(sdmat%err)
#endif
      end if
      deallocate(sdmat)
    end if
    ierr = 0

  end subroutine phist_DsdMat_delete


  subroutine phist_DsdMat_extract_view(sdmat_ptr, raw_ptr, lda, ierr) bind(C,name='phist_DsdMat_extract_view_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: sdmat_ptr
    type(C_PTR),        intent(out) :: raw_ptr
    integer(C_INT32_T), intent(out) :: lda
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

#ifdef F_DEBUG
    write(*,*) 'extract view of sdmat at address', transfer(sdmat_ptr,dummy)
    flush(6)
#endif
    if( c_associated(sdmat_ptr) ) then
      call c_f_pointer(sdmat_ptr, sdmat)
      raw_ptr = c_loc(sdmat%val(sdmat%imin,sdmat%jmin))
      lda = size(sdmat%val,1)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_DsdMat_extract_view

#ifdef PHIST_HIGH_PRECISION_KERNELS

  subroutine phist_DsdMat_extract_error(sdmat_ptr, raw_ptr, ierr) bind(C,name='phist_DsdMat_extract_error_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: sdmat_ptr
    type(C_PTR),        intent(out) :: raw_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

    if( c_associated(sdmat_ptr) ) then
      call c_f_pointer(sdmat_ptr, sdmat)
      raw_ptr = c_loc(sdmat%err(sdmat%imin,sdmat%jmin))
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_DsdMat_extract_error

#endif

  subroutine phist_DsdMat_get_nrows(sdmat_ptr, nrows, ierr) bind(C,name='phist_DsdMat_get_nrows_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: sdmat_ptr
    integer(C_INT),     intent(out) :: nrows
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

    if( c_associated(sdmat_ptr) ) then
      call c_f_pointer(sdmat_ptr, sdmat)
      nrows = sdmat%imax-sdmat%imin+1
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_DsdMat_get_nrows


  subroutine phist_DsdMat_get_ncols(sdmat_ptr, ncols, ierr) bind(C,name='phist_DsdMat_get_ncols_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: sdmat_ptr
    integer(C_INT),     intent(out) :: ncols
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

    if( c_associated(sdmat_ptr) ) then
      call c_f_pointer(sdmat_ptr, sdmat)
      ncols = sdmat%jmax-sdmat%jmin+1
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_DsdMat_get_ncols


  subroutine phist_DsdMat_view_block(sdmat_ptr, view_ptr, imin, imax, jmin, jmax, ierr) bind(C,name='phist_DsdMat_view_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    type(C_PTR),        intent(inout) :: view_ptr
    integer(C_INT),     value         :: imin, imax, jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat, view
#ifdef F_DEBUG
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

#ifdef F_DEBUG
    write(*,*) 'create view of sdmat at address', transfer(sdmat_ptr,dummy)
    flush(6)
#endif
    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(sdmat_ptr, sdmat)

    if( c_associated(view_ptr) ) then
      call c_f_pointer(view_ptr, view)
      if( .not. view%is_view ) then
        write(*,*) 'is not a view'
        flush(6)
        ierr = -88
        return
      end if
#ifdef F_DEBUG
      write(*,*) 'reusing view at address', transfer(view_ptr,dummy)
      flush(6)
#endif
    else
      allocate(view)
      view_ptr = c_loc(view)
#ifdef F_DEBUG
      write(*,*) 'created new view at address', transfer(view_ptr,dummy)
      flush(6)
#endif
    end if


    view%imin = sdmat%imin+imin
    view%imax = sdmat%imin+imax
    view%jmin = sdmat%jmin+jmin
    view%jmax = sdmat%jmin+jmax
    view%comm = sdmat%comm
    view%val => sdmat%val
#ifdef PHIST_HIGH_PRECISION_KERNELS
    view%err => sdmat%err
#endif
    view%is_view = .true.

    ierr = 0

  end subroutine phist_DsdMat_view_block


  subroutine phist_DsdMat_get_block(sdmat_ptr, block_ptr, imin, imax, jmin, jmax, ierr) bind(C,name='phist_DsdMat_get_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr, block_ptr
    integer(C_INT),     value         :: imin, imax, jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat, block
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) .or. .not. c_associated(block_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)
    call c_f_pointer(block_ptr,block)

    if( imax-imin .ne. block%imax-block%imin .or. &
      & jmax-jmin .ne. block%jmax-block%jmin .or. &
      & imax-imin .gt. sdmat%imax-sdmat%imin .or. &
      & jmax-jmin .gt. sdmat%jmax-sdmat%jmin      ) then
      ierr = -1
      return
    end if

    block%val(block%imin:block%imax,block%jmin:block%jmax) = &
      & sdmat%val(sdmat%imin+imin:sdmat%imin+imax,sdmat%jmin+jmin:sdmat%jmin+jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    block%err(block%imin:block%imax,block%jmin:block%jmax) = &
      & sdmat%err(sdmat%imin+imin:sdmat%imin+imax,sdmat%jmin+jmin:sdmat%jmin+jmax)
#endif
    ierr = 0

  end subroutine phist_DsdMat_get_block


  subroutine phist_DsdMat_set_block(sdmat_ptr, block_ptr, imin, imax, jmin, jmax, ierr) bind(C,name='phist_DsdMat_set_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr, block_ptr
    integer(C_INT),     value         :: imin, imax, jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat, block
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) .or. .not. c_associated(block_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)
    call c_f_pointer(block_ptr,block)

    if( imax-imin .ne. block%imax-block%imin .or. &
      & jmax-jmin .ne. block%jmax-block%jmin .or. &
      & imax-imin .gt. sdmat%imax-sdmat%imin .or. &
      & jmax-jmin .gt. sdmat%jmax-sdmat%jmin      ) then
      ierr = -1
      return
    end if

    sdmat%val(sdmat%imin+imin:sdmat%imin+imax,sdmat%jmin+jmin:sdmat%jmin+jmax) = &
      & block%val(block%imin:block%imax,block%jmin:block%jmax)
#ifdef PHIST_HIGH_PRECISION_KERNELS
    sdmat%err(sdmat%imin+imin:sdmat%imin+imax,sdmat%jmin+jmin:sdmat%jmin+jmax) = &
      & block%err(block%imin:block%imax,block%jmin:block%jmax)
#endif
    ierr = 0

  end subroutine phist_DsdMat_set_block


  subroutine phist_DsdMat_put_value(sdmat_ptr, val, ierr) bind(C,name='phist_DsdMat_put_value_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    real(C_DOUBLE),     value         :: val
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)

    sdmat%val(sdmat%imin:sdmat%imax,sdmat%jmin:sdmat%jmax) = val
#ifdef PHIST_HIGH_PRECISION_KERNELS
    sdmat%err(sdmat%imin:sdmat%imax,sdmat%jmin:sdmat%jmax) = 0._8
#endif
    ierr = 0

  end subroutine phist_DsdMat_put_value


  subroutine phist_DsdMat_random(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_random_f')
    use, intrinsic :: iso_c_binding
    use mpi_f08
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat
    real(kind=8), allocatable :: buff(:,:)
    integer :: m, n
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)

    m = sdmat%imax-sdmat%imin+1
    n = sdmat%jmax-sdmat%jmin+1
    allocate(buff(m,n))
    call phist_Drandom_number(size(buff), buff(1,1))

    sdmat%val(sdmat%imin:sdmat%imax,sdmat%jmin:sdmat%jmax) = buff
#ifdef PHIST_HIGH_PRECISION_KERNELS
    sdmat%err(sdmat%imin:sdmat%imax,sdmat%jmin:sdmat%jmax) = 0._8
#endif

    ierr = 0
  end subroutine phist_DsdMat_random


  subroutine phist_DsdMat_identity(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_identity_f')
    use, intrinsic :: iso_c_binding
    use mpi_f08
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat
    integer :: i, j
    !--------------------------------------------------------------------------------
    
    ierr=0

    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)

    do i = sdmat%imin, sdmat%imax, 1
      do j = sdmat%jmin, sdmat%jmax, 1
        if( i-sdmat%imin .eq. j-sdmat%jmin ) then
          sdmat%val(i,j) = 1._8
        else
          sdmat%val(i,j) = 0._8
        end if
#ifdef PHIST_HIGH_PRECISION_KERNELS
        sdmat%err(i,j) = 0._8
#endif
      end do
    end do

  end subroutine phist_DsdMat_identity


  subroutine phist_DsdMat_print(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_print_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat
    integer :: i, j
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)

    do i = sdmat%imin, sdmat%imax
      do j = sdmat%jmin, sdmat%jmax
#ifdef PHIST_HIGH_PRECISION_KERNELS
        write(*,'(G25.16,A1,G25.16)',advance='no') sdmat%val(i,j),'+',sdmat%err(i,j)
#else
        write(*,'(G16.8,A1,G25.16)',advance='no') sdmat%val(i,j)
#endif
      end do
      write(*,*)
    end do
    flush(6)
    ierr = 0

  end subroutine phist_DsdMat_print


  subroutine phist_DsdMat_add_sdMat(alpha, A_ptr, beta, B_ptr, ierr) bind(C,name='phist_DsdMat_add_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: A_ptr, B_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: A, B
    !--------------------------------------------------------------------------------

    if( .not. c_associated(A_ptr) .or. .not. c_associated(B_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(A_ptr, A)
    call c_f_pointer(B_ptr, B)

    if( A%imax-A%imin .ne. B%imax-B%imin .or. &
      & A%jmax-A%jmin .ne. B%jmax-B%jmin      ) then
      ierr = -1
      return
    end if

    call sdmat_add_sdmat(alpha, A, beta, B)

    ierr = 0

  end subroutine phist_DsdMat_add_sdMat


  subroutine phist_DsdMatT_add_sdMat(alpha, A_ptr, beta, B_ptr, ierr) bind(C,name='phist_DsdMatT_add_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: A_ptr, B_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: A, B
    !--------------------------------------------------------------------------------

    if( .not. c_associated(A_ptr) .or. .not. c_associated(B_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(A_ptr, A)
    call c_f_pointer(B_ptr, B)

    if( A%imax-A%imin .ne. B%jmax-B%jmin .or. &
      & A%jmax-A%jmin .ne. B%imax-B%imin      ) then
      ierr = -1
      return
    end if

    call sdmatT_add_sdmat(alpha, A, beta, B)

    ierr = 0

  end subroutine phist_DsdMatT_add_sdMat


  subroutine phist_DsdMat_times_sdMat(alpha, A_ptr, B_ptr, beta, M_ptr, ierr) bind(C,name='phist_DsdMat_times_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: A_ptr, B_ptr, M_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: A, B, M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(A_ptr) .or. &
      & .not. c_associated(B_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(A_ptr, A)
    call c_f_pointer(B_ptr, B)
    call c_f_pointer(M_ptr, M)

    if( A%jmax-A%jmin .ne. B%imax-B%imin .or. &
      & A%imax-A%imin .ne. M%imax-M%imin .or. &
      & B%jmax-B%jmin .ne. M%jmax-M%jmin      ) then
      ierr = -1
      return
    end if

    call sdmat_times_sdmat('N','N',alpha, A, B, beta, M, ierr)

  end subroutine phist_DsdMat_times_sdMat


  subroutine phist_DsdMatT_times_sdMat(alpha, A_ptr, B_ptr, beta, M_ptr, ierr) bind(C,name='phist_DsdMatT_times_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: A_ptr, B_ptr, M_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: A, B, M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(A_ptr) .or. &
      & .not. c_associated(B_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(A_ptr, A)
    call c_f_pointer(B_ptr, B)
    call c_f_pointer(M_ptr, M)

    if( A%imax-A%imin .ne. B%imax-B%imin .or. &
      & A%jmax-A%jmin .ne. M%imax-M%imin .or. &
      & B%jmax-B%jmin .ne. M%jmax-M%jmin      ) then
      ierr = -1
      return
    end if

    call sdmat_times_sdmat('T','N',alpha, A, B, beta, M, ierr)

  end subroutine phist_DsdMatT_times_sdMat


  subroutine phist_DsdMat_times_sdMatT(alpha, A_ptr, B_ptr, beta, M_ptr, ierr) bind(C,name='phist_DsdMat_times_sdMatT_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: A_ptr, B_ptr, M_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: A, B, M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(A_ptr) .or. &
      & .not. c_associated(B_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(A_ptr, A)
    call c_f_pointer(B_ptr, B)
    call c_f_pointer(M_ptr, M)

    if( A%imax-A%imin .ne. M%imax-M%imin .or. &
      & A%jmax-A%jmin .ne. B%jmax-B%jmin .or. &
      & B%imax-B%imin .ne. M%jmax-M%jmin      ) then
      ierr = -1
      return
    end if

    call sdmat_times_sdmat('N','T',alpha, A, B, beta, M, ierr)

  end subroutine phist_DsdMat_times_sdMatT

end module sdmat_module
