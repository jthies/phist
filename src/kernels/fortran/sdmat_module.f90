module sdmat_module
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
  !public :: phist_DsdMat_add_sdMat
  public :: sdmat_add_sdmat
  !public :: phist_DsdMat_times_sdMat
  !public :: phist_DsdMatT_times_sdMat
  public :: sdmat_times_sdmat


  !==================================================================================
  !> sdmat with row-wise layout
  type SDMat_t
    !--------------------------------------------------------------------------------
    integer     :: imin, imax
    integer     :: jmin, jmax
    integer     :: comm
    real(kind=8), contiguous, pointer :: val(:,:) => null()
    logical     :: is_view
    !--------------------------------------------------------------------------------
  end type SDMat_t

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

    B%val(B%imin:B%imax,B%jmin:B%jmax) = &
      & alpha * A%val(A%imin:A%imax,A%jmin:A%jmax)  + &
      & beta  * B%val(B%imin:B%imax,B%jmin:B%jmax)

  end subroutine sdmat_add_sdmat


  !==================================================================================
  ! gemm operation for SDMat_t
  subroutine sdmat_times_sdmat(transA,transB,alpha,A,B,beta,C)
    !--------------------------------------------------------------------------------
    character(len=1),intent(in)   :: transA, transB
    real(kind=8),   intent(in)    :: alpha
    type(SDMat_t),  intent(in)    :: A
    type(SDMat_t),  intent(in)    :: B
    real(kind=8),   intent(in)    :: beta
    type(SDMat_t),  intent(inout) :: C
    !--------------------------------------------------------------------------------
    integer :: m, n, k, lda, ldb, ldc
    !--------------------------------------------------------------------------------

    ! determine data layout
    m = C%imax-C%imin+1
    n = C%jmax-C%jmin+1
    if( transA .eq. 'N' ) then
      k = A%jmax-A%jmin+1
    else if( transA .eq. 'T' ) then
      k = A%imax-A%imin+1
    else
      call exit(1)
    end if
    lda = size(A%val,1)
    ldb = size(B%val,1)
    ldc = size(C%val,1)

    ! just call the BLAS
    call dgemm(transA,transB,m,n,k,alpha,A%val(A%imin,A%jmin),lda, &
      &                                  B%val(B%imin,B%jmin),ldb, &
      &                            beta, C%val(C%imin,C%jmin),ldc  )

    !--------------------------------------------------------------------------------
  end subroutine sdmat_times_sdmat


  !==================================================================================
  ! wrapper routines 

  subroutine phist_DsdMat_create(sdmat_ptr, nrows, ncols, comm_ptr, ierr) bind(C,name='phist_DsdMat_create_f')
    use, intrinsic :: iso_c_binding
    use mpi
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: sdmat_ptr
    integer(C_INT),     value       :: nrows, ncols
    type(C_PTR),        value       :: comm_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    integer, pointer :: comm
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
      sdmat%comm = MPI_COMM_NULL
    end if
#ifdef TESTING
    write(*,*) 'creating new sdmat with dimensions:', nrows, ncols, 'address', c_loc(sdmat)
    flush(6)
#endif
    allocate(sdmat%val(nrows,ncols))
    sdmat_ptr = c_loc(sdmat)
    ierr = 0

  end subroutine phist_DsdMat_create


  subroutine phist_DsdMat_delete(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_delete_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: sdmat_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(sdmat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

#ifdef TESTING
    write(*,*) 'deleting sdmat at address', sdmat_ptr
    flush(6)
#endif
    if( c_associated(sdmat_ptr) ) then
      call c_f_pointer(sdmat_ptr, sdmat)
      if( .not. sdmat%is_view) then
        deallocate(sdmat%val)
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
    !--------------------------------------------------------------------------------

#ifdef TESTING
    write(*,*) 'extract view of sdmat at address', sdmat_ptr
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
    !--------------------------------------------------------------------------------

#ifdef TESTING
    write(*,*) 'create view of sdmat at address', sdmat_ptr
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
#ifdef TESTING
      write(*,*) 'reusing view at address', view_ptr
      flush(6)
#endif
    else
      allocate(view)
      view_ptr = c_loc(view)
#ifdef TESTING
      write(*,*) 'created new view at address', view_ptr
      flush(6)
#endif
    end if


    view%imin = sdmat%imin+imin
    view%imax = sdmat%imin+imax
    view%jmin = sdmat%jmin+jmin
    view%jmax = sdmat%jmin+jmax
    view%comm = sdmat%comm
    view%val => sdmat%val
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

    block%val(block%imin:block%imax,block%jmin:block%jmax) = &
      & sdmat%val(sdmat%imin+imin:sdmat%imin+imax,sdmat%jmin+jmin:sdmat%jmin+jmax)
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

    sdmat%val(sdmat%imin+imin:sdmat%imin+imax,sdmat%jmin+jmin:sdmat%jmin+jmax) = &
      & block%val(block%imin:block%imax,block%jmin:block%jmax)
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
    ierr = 0

  end subroutine phist_DsdMat_put_value


  subroutine phist_DsdMat_random(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_random_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)

    call random_number(sdmat%val(sdmat%imin:sdmat%imax,sdmat%jmin:sdmat%jmax))
    ierr = 0

  end subroutine phist_DsdMat_random


  subroutine phist_DsdMat_print(sdmat_ptr, ierr) bind(C,name='phist_DsdMat_print_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: sdmat_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(SDMat_t), pointer :: sdmat
    !--------------------------------------------------------------------------------

    if( .not. c_associated(sdmat_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(sdmat_ptr, sdmat)

    write(*,*) sdmat%val(sdmat%imin:sdmat%imax,sdmat%jmin:sdmat%jmax)
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

    call sdmat_add_sdmat(alpha, A, beta, B)

    ierr = 0

  end subroutine phist_DsdMat_add_sdMat


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

    call sdmat_times_sdmat('N','N',alpha, A, B, beta, M)

    ierr = 0

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

    call sdmat_times_sdmat('T','N',alpha, A, B, beta, M)

    ierr = 0

  end subroutine phist_DsdMatT_times_sdMat



end module sdmat_module
