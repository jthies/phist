module crsmat_module
  use map_module, only: Map_t, map_setup
  use mvec_module, only: MVec_t, mvec_scale
  implicit none
  private

  public :: CrsMat_t
  !public :: phist_DcrsMat_read_mm
  !public :: phist_DcrsMat_delete
  !public :: phist_DcrsMat_times_mvec
  public :: crsmat_times_mvec
  !public :: phist_DcrsMat_get_map

  !==================================================================================
  !> simple crs matrix (not distributed!)
  type CrsMat_t
    !--------------------------------------------------------------------------------
    integer                   :: nRows                             !< number of rows
    integer                   :: nCols                             !< number of columns
    integer                   :: nEntries                          !< number of non-zero entries
    integer,      allocatable :: row_offset(:)                     !< indices of rows in col_idx and val
    integer,      allocatable :: col_idx(:)                        !< column indices
    real(kind=8), allocatable :: val(:)                            !< matrix entries
    type(Map_t)               :: row_map
    !--------------------------------------------------------------------------------
  end type CrsMat_t


contains


  !==================================================================================
  !> multiply crsmat with mvec
  subroutine crsmat_times_mvec(alpha, A, x, beta, y)
    !--------------------------------------------------------------------------------
    real(kind=8),   intent(in)    :: alpha
    type(CrsMat_t), intent(in)    :: A
    type(MVec_t),   intent(in)    :: x
    real(kind=8),   intent(in)    :: beta
    type(MVec_t),   intent(inout) :: y
    !--------------------------------------------------------------------------------
    integer :: nvec, ldx, ldy
    logical :: strided_x, strided_y
    !--------------------------------------------------------------------------------

    ! if alpha == 0, only scale y
    if( alpha .eq. 0 .or. A%nEntries .eq. 0 ) then
      call mvec_scale(y,beta)
      return
    end if

    ! determin data layout
    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided_x = .false.
    else
      strided_x = .true.
    end if

    if( .not. y%is_view .or. &
      & ( y%jmin .eq. lbound(y%val,1) .and. &
      &   y%jmax .eq. ubound(y%val,1)       ) ) then
      strided_y = .false.
    else
      strided_y = .true.
    end if
    nvec = x%jmax-x%jmin+1
    ldx = size(x%val,1)
    ldy = size(y%val,1)

    ! try to use NT stores if possible
    if( beta .eq. 0 ) then
      if( nvec .eq. 1 ) then
        if( .not. strided_x ) then
          call dspmvm_NT_1(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, x%val, y%val)
          return
        end if
      else if( nvec .eq. 2 ) then
        if( strided_x ) then
          call dspmvm_NT_strided_2(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            &                      x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call dspmvm_NT_2(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            &              x%val, y%val(y%jmin,1), ldy)
        end if
        return
      else if( nvec .eq. 4 ) then
        if( strided_x ) then
          call dspmvm_NT_strided_4(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            &                      x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call dspmvm_NT_4(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            &              x%val, y%val(y%jmin,1), ldy)
        end if
        return
      else if( nvec .eq. 8 ) then
        if( strided_x ) then
          call dspmvm_NT_strided_8(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            &                      x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call dspmvm_NT_8(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            &              x%val, y%val(y%jmin,1), ldy)
        end if
        return
      end if
    end if

    if( nvec .eq. 1 ) then
      if( strided_x .or. strided_y ) then
        call dspmvm_strided_1(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &                   x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call dspmvm_1(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &           x%val, beta, y%val)
      end if
      return
    else if( nvec .eq. 2 ) then
      if( strided_x .or. strided_y ) then
        call dspmvm_strided_2(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &                   x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call dspmvm_2(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &           x%val, beta, y%val)
      end if
      return
    else if( nvec .eq. 4 ) then
      if( strided_x .or. strided_y ) then
        call dspmvm_strided_4(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &                   x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call dspmvm_4(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &           x%val, beta, y%val)
      end if
      return
    else if( nvec .eq. 8 ) then
      if( strided_x .or. strided_y ) then
        call dspmvm_strided_8(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &                   x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call dspmvm_8(A%nrows, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          &           x%val, beta, y%val)
      end if
      return
    end if

    call dspmvm_generic(A%nrows, nvec, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
      &                 x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
    !--------------------------------------------------------------------------------
  end subroutine crsmat_times_mvec


  !==================================================================================
  !> read MatrixMarket file
  subroutine phist_DcrsMat_read_mm(A_ptr, filename_len, filename_ptr, ierr) bind(C,name='phist_DcrsMat_read_mm_f')
    use, intrinsic :: iso_c_binding
    use m_mrgrnk
    use env_module, only: newunit
    use mpi
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: A_ptr
    integer(C_INT),     value       :: filename_len
    character(C_CHAR),  intent(in)  :: filename_ptr(filename_len)
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    character(len=filename_len) :: filename
    !--------------------------------------------------------------------------------
    integer :: funit
    character(len=100) :: line
    integer, allocatable :: idx(:,:)
    real(kind=8),    allocatable :: val(:)
    integer :: i, j, off, n
    !--------------------------------------------------------------------------------

    do i = 1, filename_len
      filename(i:i) = filename_ptr(i)
    end do

    ! open the file
    write(*,*) 'reading file:', filename
    open(unit   = newunit(funit), file    = filename, &
      &  action = 'read',         status  = 'old',    &
      &  iostat = ierr)
    if( ierr .ne. 0 ) return

    ! read first line
    read(funit,'(A)') line
    write(*,*) line
    if( trim(line) .ne. '%%MatrixMarket matrix coordinate real general' ) then
      write(*,*) 'unsupported format'
      ierr = -99
      return
    end if

    ! read second line (comment)
    read(funit,*) line
    if( line(1:1) .ne. '%' ) then
      call abort
    end if

    allocate(A)

    ! now read the dimensions
    read(funit,*) A%nRows, A%nCols, A%nEntries
    write(*,*) 'CrsMat:', A%nRows, A%nCols, A%nEntries

    call map_setup(A%row_map, MPI_COMM_WORLD, int(A%nRows,kind=8), ierr)
    if( ierr .ne. 0 ) return

    ! allocate temporary buffers
    allocate(idx(A%nEntries,2))
    allocate(val(A%nEntries))

    ! read data
    do i = 1, A%nEntries
      read(funit,*) idx(i,1), idx(i,2), val(i)
    end do

    ! close the file
    close(funit)


    ! allocate crs matrix
    allocate(A%row_offset(A%nRows+1))
    allocate(A%col_idx(A%nEntries))
    allocate(A%val(A%nEntries))

    ! count number of entries per row
    A%row_offset = 0
    do i = 1, A%nEntries
      A%row_offset( idx(i,1)+1 ) = A%row_offset( idx(i,1)+1 )  +  1
    end do
    A%row_offset(1) = 1
    do i = 1, A%nRows
      A%row_offset( i+1 ) = A%row_offset( i ) + A%row_offset( i+1 )
    end do

    ! now put the entries into the corresponding rows
    A%col_idx = 0
    A%val = 0.
    do i = 1, A%nEntries
      A%col_idx( A%row_offset(idx(i,1)) ) = idx(i,2)
      A%val    ( A%row_offset(idx(i,1)) ) = val(i)
      A%row_offset( idx(i,1) ) = A%row_offset( idx(i,1) )  +  1
    end do
    ! set row_offset
    do i = A%nRows, 1, -1
      A%row_offset( i+1 ) = A%row_offset( i )
    end do
    A%row_offset( 1 ) = 1

    ! we still need to sort the entries in each row by column
    do i = 1, A%nRows
      off = A%row_offset(i)
      n = A%row_offset(i+1) - off
      if( n .gt. 1 ) then
        ! reuse idx for sorting
        idx(1:n,1) = A%col_idx( off:off+n-1 )
        val(1:n)   = A%val( off:off+n-1 )

        ! determine order of columns
        call mrgrnk(idx(1:n,1), idx(1:n,2))

        ! copy back sorted arrays
        do j = 1, n
          A%col_idx( off+idx(j,2)-1 ) = idx(j,1)
          A%val( off+idx(j,2)-1 ) = val(j)
        end do
      end if
    end do

    write(*,*) 'created new crsMat with dimensions', A%nRows, A%nCols, A%nEntries, 'address', c_loc(A)
    A_ptr = c_loc(A)

    ierr = 0

    !--------------------------------------------------------------------------------
  end subroutine phist_DcrsMat_read_mm


  !==================================================================================
  subroutine phist_DcrsMat_delete(A_ptr, ierr) bind(C,name='phist_DcrsMat_delete_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),    value       :: A_ptr
    integer(C_INT), intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    !--------------------------------------------------------------------------------

    write(*,*) 'deleting crsMat at address', A_ptr
    if( c_associated(A_ptr) ) then
      call c_f_pointer(A_ptr,A)
      deallocate(A)
    end if
    ierr = 0

    !--------------------------------------------------------------------------------
  end subroutine phist_DcrsMat_delete


  !==================================================================================
  subroutine phist_DcrsMat_get_map(A_ptr, map_ptr, ierr) bind(C,name='phist_DcrsMat_get_map_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),    value       :: A_ptr
    type(C_PTR),    intent(out) :: map_ptr
    integer(C_INT), intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    !--------------------------------------------------------------------------------

    if( c_associated(A_ptr) ) then
      call c_f_pointer(A_ptr,A)
      map_ptr = c_loc(A%row_map)
    else
      ierr = -88
    end if

    !--------------------------------------------------------------------------------
  end subroutine phist_DcrsMat_get_map


  !==================================================================================
  subroutine phist_DcrsMat_times_mvec(alpha, A_ptr, x_ptr, beta, y_ptr, ierr) bind(C,name='phist_DcrsMat_times_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    real(C_DOUBLE),   value         :: alpha, beta
    type(C_PTR),      value         :: A_ptr, x_ptr, y_ptr
    integer(C_INT),   intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    type(MVec_t), pointer :: x, y
    !--------------------------------------------------------------------------------

    if( .not. c_associated(A_ptr) .or. &
      & .not. c_associated(x_ptr) .or. &
      & .not. c_associated(y_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(A_ptr,A)
    call c_f_pointer(x_ptr,x)
    call c_f_pointer(y_ptr,y)

    call crsmat_times_mvec(alpha,A,x,beta,y)

    ierr = 0
  end subroutine phist_DcrsMat_times_mvec


end module crsmat_module

