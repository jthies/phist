module mvec_module
  use map_module, only: Map_t
  use sdmat_module, only: SDMat_t
  implicit none
  private

  public :: MVec_t
  !public :: phist_Dmvec_create
  !public :: phist_Dmvec_delete
  !public :: phist_Dmvec_extract_view
  !public :: phist_Dmvec_my_length
  !public :: phist_Dmvec_get_map
  !public :: phist_Dmvec_num_vectors
  !public :: phist_Dmvec_view_block
  !public :: phist_Dmvec_get_block
  !public :: phist_Dmvec_set_block
  !public :: phist_Dmvec_gather_mvecs
  !public :: phist_Dmvec_scatter_mvecs
  public :: mvec_gather_mvecs
  public :: mvec_scatter_mvecs
  !public :: phist_Dmvec_put_value
  !public :: phist_Dmvec_random
  !public :: phist_Dmvec_print
  !public :: phist_Dmvec_norm2
  public :: mvec_norm2
  !public :: phist_Dmvec_scale
  !public :: phist_Dmvec_vscale
  public :: mvec_scale
  public :: mvec_vscale
  !public :: phist_Dmvec_add_mvec
  !public :: phist_Dmvec_vadd_mvec
  public :: mvec_add_mvec
  public :: mvec_vadd_mvec
  !public :: phist_Dmvec_dot_mvec
  public :: mvec_dot_mvec
  !public :: phist_Dmvec_times_sdMat
  !public :: phist_Dmvec_times_sdMat_inplace
  public :: mvec_times_sdmat
  public :: mvec_times_sdmat_inplace
  !public :: phist_DmvecT_times_mvec
  public :: mvecT_times_mvec
  !public :: phist_Dmvec_QR
  public :: mvec_QR


  !==================================================================================
  !> mvec with row-wise layout
  type MVec_t
    !--------------------------------------------------------------------------------
    integer     :: jmin, jmax
    type(Map_t) :: map
    real(kind=8), contiguous, pointer :: val(:,:) => null()
    logical     :: is_view
    !--------------------------------------------------------------------------------
  end type MVec_t

contains

  !==================================================================================
  !> copy data from a list of other mvecs into one mvec
  subroutine mvec_gather_mvecs(mvec, block_list)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: mvec
    type(MVec_t), intent(in)    :: block_list(1:)
    !--------------------------------------------------------------------------------
    integer :: nvec, nblocks, nrows, ldm
    logical :: single_vector_gather
    integer :: i, j, nvec_i, off, jmin_i, jmax_i
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    nvec = mvec%jmax-mvec%jmin+1
    nblocks = size(block_list)
    nrows = size(mvec%val,2)
    ldm = size(mvec%val,1)
    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if


    ! fast versions for fixed nvec and gathering single vectors
    single_vector_gather = .true.
    if( nvec .ne. nblocks ) single_vector_gather = .false.
    do i = 1, nblocks
      if( block_list(i)%jmin .ne. block_list(i)%jmax ) then
        single_vector_gather = .false.
      end if
    end do

    if( single_vector_gather .and. .not. strided ) then
      if( nvec .eq. 1 ) then
        call dgather_1(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1) )
        return
      else if( nvec .eq. 2 ) then
        call dgather_2(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(2)%val(block_list(2)%jmin,1), size(block_list(2)%val,1) )
        return
      else if( nvec .eq. 4 ) then
        call dgather_4(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(2)%val(block_list(2)%jmin,1), size(block_list(2)%val,1), &
          &            block_list(3)%val(block_list(3)%jmin,1), size(block_list(3)%val,1), &
          &            block_list(4)%val(block_list(4)%jmin,1), size(block_list(4)%val,1) )
        return
      end if
    end if

!$omp parallel do private(off) schedule(static)
    do j = 1, nrows
      off = mvec%jmin
      do i = 1, nblocks
        jmin_i = block_list(i)%jmin
        jmax_i = block_list(i)%jmax
        nvec_i = jmax_i-jmin_i+1
        mvec%val(off:off+nvec_i-1,j) = block_list(i)%val(jmin_i:jmax_i,j)
        off = off +nvec_i
      end do
    end do

    !--------------------------------------------------------------------------------
  end subroutine mvec_gather_mvecs


  ! currently unused
  !==================================================================================
  !> copy data from a list of other mvecs into one mvec
  subroutine mvec_scatter_mvecs(mvec, block_list)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: mvec
    type(MVec_t), intent(in)    :: block_list(1:)
    !--------------------------------------------------------------------------------
    integer :: nvec, nblocks, nrows, ldm
    logical :: single_vector_scatter
    integer :: i, j, nvec_i, off, jmin_i, jmax_i
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    nvec = mvec%jmax-mvec%jmin+1
    nblocks = size(block_list)
    nrows = size(mvec%val,2)
    ldm = size(mvec%val,1)
    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if


    ! fast versions for fixed nvec and scattering single vectors
    single_vector_scatter = .true.
    if( nvec .ne. nblocks ) single_vector_scatter = .false.
    do i = 1, nblocks
      if( block_list(i)%jmin .ne. block_list(i)%jmax ) then
        single_vector_scatter = .false.
      end if
    end do

    if( single_vector_scatter .and. .not. strided ) then
      if( nvec .eq. 1 ) then
        call dscatter_1(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1) )
        return
      else if( nvec .eq. 2 ) then
        call dscatter_2(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(1)%val(block_list(2)%jmin,1), size(block_list(2)%val,1) )
        return
      else if( nvec .eq. 4 ) then
        call dscatter_4(nrows, mvec%val(mvec%jmin,1), &
          &            block_list(1)%val(block_list(1)%jmin,1), size(block_list(1)%val,1), &
          &            block_list(2)%val(block_list(2)%jmin,1), size(block_list(2)%val,1), &
          &            block_list(3)%val(block_list(3)%jmin,1), size(block_list(3)%val,1), &
          &            block_list(4)%val(block_list(4)%jmin,1), size(block_list(4)%val,1) )
        return
      end if
    end if

!$omp parallel do private(off) schedule(static)
    do j = 1, nrows
      off = mvec%jmin
      do i = 1, nblocks
        jmin_i = block_list(i)%jmin
        jmax_i = block_list(i)%jmax
        nvec_i = jmax_i-jmin_i+1
        block_list(i)%val(jmin_i:jmax_i,j) = mvec%val(off:off+nvec_i-1,j)
        off = off +nvec_i
      end do
    end do

    !--------------------------------------------------------------------------------
  end subroutine mvec_scatter_mvecs



  !==================================================================================
  !> calculate 2-norm of the vectors in a multivector
  subroutine mvec_norm2(mvec, vnrm)
    use mpi
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(in)  :: mvec
    real(kind=8), intent(out) :: vnrm(mvec%jmin:mvec%jmax)
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, lda, ierr
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    if( .not. mvec%is_view .or. &
      & ( mvec%jmin .eq. lbound(mvec%val,1) .and. &
      &   mvec%jmax .eq. ubound(mvec%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if

    nvec = mvec%jmax-mvec%jmin+1
    nrows = size(mvec%val,2)
    lda = size(mvec%val,1)
    ! for single vectors call appropriate blas
    if( nvec .eq. 1 ) then
      call dnrm2_strided_1(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call dnrm2_strided_2(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_2(nrows, mvec%val, vnrm)
      end if
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call dnrm2_strided_4(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_4(nrows, mvec%val, vnrm)
      end if
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call dnrm2_strided_8(nrows, mvec%val(mvec%jmin,1), lda, vnrm)
      else
        call dnrm2_8(nrows, mvec%val, vnrm)
      end if
    else
      call dnrm2_general(nrows, nvec, mvec%val(mvec%jmin,1), lda, vnrm)
    end if

    vnrm = vnrm*vnrm
    call MPI_Allreduce(MPI_IN_PLACE,vnrm,nvec,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    vnrm = sqrt(vnrm)

    !--------------------------------------------------------------------------------
  end subroutine mvec_norm2


  !==================================================================================
  ! scale mvec
  subroutine mvec_scale(x,alpha)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: x
    real(kind=8), intent(in)    :: alpha
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows
    logical :: strided
    real(kind=8), allocatable :: alpha_vec(:)
    !--------------------------------------------------------------------------------

    if( alpha .eq. 1 ) return

    ! determine data layout
    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if

    nvec = x%jmax-x%jmin+1
    nrows = size(x%val,2)

    if( .not. strided ) then
      call dscal(nvec*nrows,alpha,x%val,1)
    else
      allocate(alpha_vec(nvec))
      alpha_vec = alpha
      call mvec_vscale(x,alpha_vec)
    end if
    !--------------------------------------------------------------------------------
  end subroutine mvec_scale


  !==================================================================================
  ! scale mvec
  subroutine mvec_vscale(x,alpha)
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(inout) :: x
    real(kind=8), intent(in)    :: alpha(*)
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, lda
    logical :: strided
    !--------------------------------------------------------------------------------

    ! determine data layout
    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided = .false.
    else
      strided = .true.
    end if

    nvec = x%jmax-x%jmin+1
    nrows = size(x%val,2)
    lda = size(x%val,1)

    !call dlascl2(nvec,nrows,alpha(1),x%val(x%jmin,1),lda)

    if( nvec .eq. 1 ) then
      if( strided ) then
        call dscal_strided_1(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_1(nrows,alpha,x%val)
      end if
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call dscal_strided_2(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_2(nrows,alpha,x%val)
      end if
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call dscal_strided_4(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_4(nrows,alpha,x%val)
      end if
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call dscal_strided_8(nrows,alpha,x%val(x%jmin,1),lda)
      else
        call dscal_8(nrows,alpha,x%val)
      end if
    else
      call dscal_general(nrows, nvec, alpha, x%val(x%jmin,1), lda)
    end if
    !--------------------------------------------------------------------------------
  end subroutine mvec_vscale


  !==================================================================================
  ! daxpy routine for mvecs, delegate to mvec_vadd_mvec
  subroutine mvec_add_mvec(alpha,x,beta,y)
    !--------------------------------------------------------------------------------
    real(kind=8), intent(in)    :: alpha
    type(MVec_t), intent(in)    :: x
    real(kind=8), intent(in)    :: beta
    type(MVec_t), intent(inout) :: y
    !--------------------------------------------------------------------------------
    real(kind=8), allocatable :: alpha_vec(:)
    !--------------------------------------------------------------------------------

    if( alpha .eq. 0 ) then
      call mvec_scale(y,beta)
    else
      allocate(alpha_vec(y%jmin:y%jmax))
      alpha_vec = alpha
      call mvec_vadd_mvec(alpha_vec,x,beta,y)
    end if

    !--------------------------------------------------------------------------------
  end subroutine mvec_add_mvec


  !==================================================================================
  ! daxpy routine for mvecs
  subroutine mvec_vadd_mvec(alpha,x,beta,y)
    !--------------------------------------------------------------------------------
    real(kind=8), intent(in)    :: alpha(*)
    type(MVec_t), intent(in)    :: x
    real(kind=8), intent(in)    :: beta
    type(MVec_t), intent(inout) :: y
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, ldx, ldy
    logical :: strided_x, strided_y, strided
    logical :: only_scale, only_copy, one_alpha, y_aligned
    integer :: i
    !--------------------------------------------------------------------------------

    ! first gather some data to decide what we want
    nvec = y%jmax-y%jmin+1
    nrows = size(y%val,2)
    ldy = size(y%val,1)

    if( .not. y%is_view .or. &
      & ( y%jmin .eq. lbound(y%val,1) .and. &
      &   y%jmax .eq. ubound(y%val,1)       ) ) then
      strided_y = .false.
    else
      strided_y = .true.
    end if


    only_scale = .true.
    do i = 1, nvec
      if( alpha(i) .ne. 0 ) only_scale = .false.
    end do
    if( only_scale ) then
      call mvec_scale(y,beta)
      return
    end if

    only_copy = .true.
    if( beta .ne. 0 ) only_copy = .false.
    do i = 1, nvec
      if( alpha(i) .ne. 1 ) only_copy = .false.
    end do

    ldx = size(x%val,1)

    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided_x = .false.
    else
      strided_x = .true.
    end if

    one_alpha = .true.
    do i = 1, nvec-1, 1
      if( alpha(i) .ne. alpha(i+1) ) one_alpha = .false.
    end do

    y_aligned = .true.
    if( mod(loc(y%val(y%jmin,1)),16) .ne. 0 .or. mod(ldy,2) .ne. 0 ) then
      y_aligned = .false.
    end if

    strided = strided_x .or. strided_y

    if( only_copy ) then
      if( .not. strided ) then
        call dcopy_1(nvec*nrows, x%val, y%val)
      else
        call dcopy_general(nvec,nrows,x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
      end if
      return
    end if


    if( beta .eq. 0 ) then
      if( nvec .eq. 2 .and. y_aligned ) then
        if( strided ) then
          call daxpy_NT_strided_2(nrows, alpha(1), x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call daxpy_NT_2(nrows, alpha(1), x%val(x%jmin,1), y%val(y%jmin,1))
        end if
        return
      else if( nvec .eq. 4 .and. y_aligned ) then
        if( strided ) then
          call daxpy_NT_strided_4(nrows, alpha(1), x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call daxpy_NT_4(nrows, alpha(1), x%val(x%jmin,1), y%val(y%jmin,1))
        end if
        return
      else if( nvec .eq. 8 .and. y_aligned ) then
        if( strided ) then
          call daxpy_NT_strided_8(nrows, alpha(1), x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy)
        else
          call daxpy_NT_8(nrows, alpha(1), x%val(x%jmin,1), y%val(y%jmin,1))
        end if
        return
      end if
    end if

    if( nvec .eq. 1 ) then
      if( strided ) then
        call daxpby_strided_1(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_1(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call daxpby_strided_2(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_2(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call daxpby_strided_4(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_4(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call daxpby_strided_8(nrows, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      else
        call daxpby_8(nrows, alpha(1), x%val(x%jmin,1), beta, y%val(y%jmin,1))
      end if
      return
    end if

    call daxpby_generic(nrows, nvec, alpha(1), x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
    !--------------------------------------------------------------------------------
  end subroutine mvec_vadd_mvec


  !==================================================================================
  ! dot product for mvecs
  subroutine mvec_dot_mvec(x,y,dot)
    use mpi
    !--------------------------------------------------------------------------------
    type(MVec_t), intent(in)  :: x, y
    real(kind=8), intent(out) :: dot(x%jmin:x%jmax)
    !--------------------------------------------------------------------------------
    integer :: nvec, nrows, ldx, ldy, ierr
    logical :: strided_x, strided_y, strided
    !--------------------------------------------------------------------------------

    ! determine data layout
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

    strided = strided_x .or. strided_y

    nvec = x%jmax-x%jmin+1
    nrows = size(x%val,2)
    ldx = size(x%val,1)
    ldy = size(y%val,1)
    ! for single vectors call appropriate blas
    if( nvec .eq. 1 ) then
      call ddot_strided_1(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
    else if( nvec .eq. 2 ) then
      if( strided ) then
        call ddot_strided_2(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_2(nrows, x%val, y%val, dot)
      end if
    else if( nvec .eq. 4 ) then
      if( strided ) then
        call ddot_strided_4(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_4(nrows, x%val, y%val, dot)
      end if
    else if( nvec .eq. 8 ) then
      if( strided ) then
        call ddot_strided_8(nrows, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
      else
        call ddot_8(nrows, x%val, y%val, dot)
      end if
    else
      call ddot_general(nrows, nvec, x%val(x%jmin,1), ldx, y%val(y%jmin,1), ldy, dot)
    end if

    call MPI_Allreduce(MPI_IN_PLACE,dot,nvec,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    !--------------------------------------------------------------------------------
  end subroutine mvec_dot_mvec


  !==================================================================================
  ! special gemm routine for mvec times sdmat
  subroutine mvec_times_sdmat(alpha,v,M,beta,w)
    !--------------------------------------------------------------------------------
    real(kind=8),  intent(in)    :: alpha, beta
    type(MVec_t),  intent(in)    :: v
    type(SDMat_t), intent(in)    :: M
    type(Mvec_t),  intent(inout) :: w
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw
    logical :: strided_v, strided_w
    real(kind=8), allocatable :: Mtmp(:,:)
    !--------------------------------------------------------------------------------

    ! check if we only need to scale
    if( alpha .eq. 0 ) then
      call mvec_scale(w,beta)
      return
    end if

    ! determine data layout
    nrows = size(w%val,2)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if

    allocate(Mtmp(nvecw,nvecv))
    Mtmp = transpose(M%val(M%imin:M%imax,M%jmin:M%jmax))
    ! recognize small block mvecs
    if( .not. strided_w ) then
      if( nvecw .eq. 1 ) then
        if( strided_v ) then
          call dgemm_sB_strided_1_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_1_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecw .eq. 2 ) then
        if( strided_v ) then
          call dgemm_sB_strided_2_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_2_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecw .eq. 4 ) then
        if( strided_v ) then
          call dgemm_sB_strided_4_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_4_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecw .eq. 8 ) then
        if( strided_v ) then
          call dgemm_sB_strided_8_k(nrows, nvecv, alpha, v%val(v%jmin,1), ldv, Mtmp, beta, w%val)
        else
          call dgemm_sB_8_k        (nrows, nvecv, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      end if
    end if
    deallocate(Mtmp)
    allocate(Mtmp(nvecv,nvecw))
    Mtmp = M%val(M%imin:M%imax,M%jmin:M%jmax)

    if( .not. strided_v ) then
      if( nvecv .eq. 1 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_1(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_1        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecv .eq. 2 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_2(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_2        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecv .eq. 4 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_4(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_4        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      else if( nvecv .eq. 8 ) then
        if( strided_w ) then
          call dgemm_sB_strided_k_8(nrows, nvecw, alpha, v%val, Mtmp, beta, w%val(w%jmin,1), ldw)
        else
          call dgemm_sB_k_8        (nrows, nvecw, alpha, v%val, Mtmp, beta, w%val)
        end if
        return
      end if
    end if


    call dgemm_sB_generic(nrows,nvecw,nvecv,alpha,v%val(v%jmin,1),ldv, Mtmp, beta, w%val(w%jmin,1),ldw)


    !--------------------------------------------------------------------------------
  end subroutine mvec_times_sdmat


  !==================================================================================
  ! special gemm routine for mvec <- mvec*sdMat
  subroutine mvec_times_sdmat_inplace(v,M)
    !--------------------------------------------------------------------------------
    type(MVec_t),  intent(in)    :: v
    type(SDMat_t), intent(in)    :: M
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv
    logical :: strided_v
    !--------------------------------------------------------------------------------

    ! determine data layout
    nrows = size(v%val,2)
    nvecv = v%jmax-v%jmin+1
    nvecw = M%jmax-M%jmin+1
    ldv = size(v%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    call dgemm_sB_generic_inplace(nrows,nvecw,nvecv,v%val(v%jmin,1),ldv, &
      &                   M%val(M%imin:M%imax,M%jmin:M%jmax))


    !--------------------------------------------------------------------------------
  end subroutine mvec_times_sdmat_inplace


  !==================================================================================
  ! special gemm routine for mvecT_times_mvec
  subroutine mvecT_times_mvec(alpha,v,w,beta,m)
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),  intent(in)    :: alpha
    type(MVec_t),  intent(in)    :: v
    type(MVec_t),  intent(in)    :: w
    real(kind=8),  intent(in)    :: beta
    type(SDMat_t), intent(inout) :: M
    !--------------------------------------------------------------------------------
    integer :: nrows, nvecv, nvecw, ldv, ldw, ierr
    logical :: strided_v, strided_w
    logical :: handled, tmp_transposed
    real(kind=8), allocatable :: tmp(:,:)
    real(kind=8), allocatable :: tmp_(:,:)
    !--------------------------------------------------------------------------------

    ! check if we only need to scale
    if( alpha .eq. 0 ) then
      M%val(M%imin:M%imax,M%jmin:M%jmax) = beta*M%val(M%imin:M%imax,M%jmin:M%jmax)
      return
    end if

    ! determine data layout
    nrows = size(w%val,2)
    nvecv = v%jmax-v%jmin+1
    nvecw = w%jmax-w%jmin+1
    ldv = size(v%val,1)
    ldw = size(w%val,1)
    if( .not. v%is_view .or. &
      & ( v%jmin .eq. lbound(v%val,1) .and. &
      &   v%jmax .eq. ubound(v%val,1)       ) ) then
      strided_v = .false.
    else
      strided_v = .true.
    end if

    if( .not. w%is_view .or. &
      & ( w%jmin .eq. lbound(w%val,1) .and. &
      &   w%jmax .eq. ubound(w%val,1)       ) ) then
      strided_w = .false.
    else
      strided_w = .true.
    end if

    handled = .false.
    tmp_transposed = .false.

    if( .not. strided_v ) then
      if( nvecv .eq. 1 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_1(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else
          call dgemm_sC_1(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      else if( nvecv .eq. 2 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_2(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else
          call dgemm_sC_2(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      else if( nvecv .eq. 4 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_4(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else
          call dgemm_sC_4(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      else if( nvecv .eq. 8 ) then
        allocate(tmp(nvecv,nvecw))
        if( strided_w ) then
          call dgemm_sC_strided_8(nrows,nvecw,v%val,w%val(w%jmin,1),ldw,tmp)
        else
          call dgemm_sC_8(nrows,nvecw,v%val,w%val,tmp)
        end if
        handled = .true.
      end if
    end if

    if( .not. handled .and. .not. strided_w ) then
      if( nvecw .eq. 1 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_w ) then
          call dgemm_sC_strided_1(nrows,nvecv,w%val,v%val(w%jmin,1),ldv,tmp)
        else
          call dgemm_sC_1(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      else if( nvecw .eq. 2 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_w ) then
          call dgemm_sC_strided_2(nrows,nvecv,w%val,v%val(w%jmin,1),ldv,tmp)
        else
          call dgemm_sC_2(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      else if( nvecw .eq. 4 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_w ) then
          call dgemm_sC_strided_4(nrows,nvecv,w%val,v%val(w%jmin,1),ldv,tmp)
        else
          call dgemm_sC_4(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      else if( nvecw .eq. 8 ) then
        allocate(tmp(nvecw,nvecv))
        if( strided_w ) then
          call dgemm_sC_strided_8(nrows,nvecv,w%val,v%val(w%jmin,1),ldv,tmp)
        else
          call dgemm_sC_8(nrows,nvecv,w%val,v%val,tmp)
        end if
        handled = .true.
        tmp_transposed = .true.
      end if
    end if


    if( .not. handled ) then
      ! generic case
      allocate(tmp(nvecv,nvecw))

      call dgemm_sC_generic(nrows,nvecv,nvecw,v%val(v%jmin,1),ldv,w%val(w%jmin,1),ldw,tmp)

    end if

    allocate(tmp_(nvecv,nvecw))
    if( tmp_transposed ) then
      tmp_ = transpose(tmp)
    else
      tmp_ = tmp
    end if
    call MPI_Allreduce(MPI_IN_PLACE,tmp_,nvecv*nvecw,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    M%val(M%imin:M%imax,M%jmin:M%jmax) = alpha*tmp_+beta*M%val(M%imin:M%imax,M%jmin:M%jmax)
    !--------------------------------------------------------------------------------
  end subroutine mvecT_times_mvec


  !==================================================================================
  ! orthogonalize v with v = QR, fill with orthogonal random vectors if not full rank
  subroutine mvec_QR(v,R,nullSpaceDim)
    !--------------------------------------------------------------------------------
    type(MVec_t),   intent(inout) :: v
    type(SDMat_t),  intent(inout) :: R
    integer,        intent(out)   :: nullSpaceDim
    !--------------------------------------------------------------------------------
    real(kind=8), parameter :: eps = 1.e-10
    integer :: i, i_, nvec, rank
    type(MVec_t) :: vi, vipn
    type(SDMat_t) :: Ripn
    real(kind=8) :: rii(1:1)
    !--------------------------------------------------------------------------------

    nvec = v%jmax-v%jmin+1

    ! setup views vi, vipn, Ripn
    vipn%jmax = v%jmax
    vipn%is_view = .true.
    vi%is_view = .true.
    vi%val=>v%val
    vipn%val=>v%val
    vi%map = v%map
    vipn%map = v%map
    Ripn%jmax = R%jmax
    Ripn%is_view = .true.
    Ripn%val => R%val
    Ripn%comm = R%comm

    ! simple modified gram schmidt... probably SLOW
    ! R = 0
    nullSpaceDim = 0
    R%val(R%imin:R%imax,R%jmin:R%jmax) = 0.
    i_ = 0
    do i = 0, nvec-1, 1
      ! create view of column i
      vi%jmin = v%jmin+i
      vi%jmax = v%jmin+i
      call mvec_norm2(vi,rii)
      R%val(R%imin+i_,R%jmin+i) = rii(1)
      if( rii(1) .lt. eps ) then
        nullSpaceDim = nullSpaceDim + 1
        rii(1) = 0
        R%val(R%imin+i_,R%jmin+i) = 0
      else
        rii(1) = 1._8/rii(1)
        call mvec_scale(vi,rii(1))
        if( i .lt. nvec-1 ) then
          vipn%jmin = v%jmin+i+1
          Ripn%imin = R%imin+i_
          Ripn%imax = R%imin+i_
          Ripn%jmin = R%jmin+i+1
          call mvecT_times_mvec(1._8,vi,vipn,0._8,Ripn)
          call mvec_times_sdmat(-1._8,vi,Ripn,1._8,vipn)
        end if
        ! copy vi to column i_, because previous vector was not linearly independent
        if( i .ne. i_ ) then
          v%val(v%jmin+i_,:) = v%val(v%jmin+i,:)
        end if
        i_ = i_ + 1
      end if
    end do

    ! try to generate orthogonal random vectors
    if( nullSpaceDim .gt. 0 ) then
      rank = nvec - nullSpaceDim
      call random_number(v%val(v%jmin+rank:v%jmax,:))

      ! reuse Ripn for temporary storage
      Ripn%val=>null()
      Ripn%is_view = .false.
      Ripn%imin = 1
      Ripn%imax = max(1,rank)
      Ripn%jmin = 1
      Ripn%jmax = nullSpaceDim
      allocate(Ripn%val(Ripn%imax,Ripn%jmax))
      Ripn%val = 0._8

      if( rank .gt. 0 ) then
        ! orthog. random vectors wrt. previous vectors
        vi%jmin = v%jmin
        vi%jmax = v%jmin+rank-1
        vipn%jmin = v%jmin+rank
        call mvecT_times_mvec(1._8,vi,vipn,0._8,Ripn)
        call mvec_times_sdmat(-1._8,vi,Ripn,1._8,vipn)
      end if

      ! orthog. random vectors wrt. each other
      Ripn%imin = 1
      Ripn%imax = 1
      do i = rank, nvec-1, 1
        ! create view of column i
        vi%jmin = v%jmin+i
        vi%jmax = v%jmin+i
        call mvec_norm2(vi,rii)
        if( rii(1) .lt. eps ) then
          write(*,*) 'error during orthogonalization'
          flush(6)
          nullSpaceDim = -1
          deallocate(Ripn%val)
          return
        end if
        rii(1) = 1._8/rii(1)
        call mvec_scale(vi,rii(1))
        if( i .lt. nvec-1 ) then
          vipn%jmin = v%jmin+i+1
          Ripn%jmin = i-rank+2
          call mvecT_times_mvec(1._8,vi,vipn,0._8,Ripn)
          call mvec_times_sdmat(-1._8,vi,Ripn,1._8,vipn)
        end if
      end do

      deallocate(Ripn%val)
    end if

    !--------------------------------------------------------------------------------
  end subroutine mvec_QR
 

  !==================================================================================
  ! wrapper routines 

  subroutine phist_Dmvec_create(mvec_ptr, map_ptr, nvec, ierr) bind(C,name='phist_Dmvec_create_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: mvec_ptr
    type(C_PTR),        value       :: map_ptr
    integer(C_INT32_T), value       :: nvec
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    type(Map_t), pointer :: map
    integer :: i
    !--------------------------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    allocate(mvec)
    mvec%is_view = .false.
    mvec%jmin = 1
    mvec%jmax = nvec
    mvec%map = map
#ifdef TESTING
    write(*,*) 'creating new mvec with dimensions:', nvec, map%nlocal(map%me), 'address', c_loc(mvec)
    flush(6)
    flush(6)
#endif
    allocate(mvec%val(nvec,map%nlocal(map%me)))
    ! that should hopefully help in cases of NUMA
!$omp parallel do schedule(static)
    do i = 1, size(mvec%val,2)
      mvec%val(:,i) = 0._8
    end do
    mvec_ptr = c_loc(mvec)
    ierr = 0

  end subroutine phist_Dmvec_create


  subroutine phist_Dmvec_delete(mvec_ptr, ierr) bind(C,name='phist_Dmvec_delete_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

#ifdef TESTING
    write(*,*) 'deleting mvec at address', mvec_ptr
    flush(6)
    flush(6)
#endif
    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      if( .not. mvec%is_view) then
        deallocate(mvec%val)
      end if
      deallocate(mvec)
    end if
    ierr = 0

  end subroutine phist_Dmvec_delete


  subroutine phist_Dmvec_extract_view(mvec_ptr, raw_ptr, lda, ierr) bind(C,name='phist_Dmvec_extract_view_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    type(C_PTR),        intent(out) :: raw_ptr
    integer(C_INT32_T), intent(out) :: lda!, stride
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

#ifdef TESTING
    !write(*,*) 'extract view of mvec at address', mvec_ptr
    !flush(6)
#endif
    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      raw_ptr = c_loc(mvec%val(mvec%jmin,1))
      lda = size(mvec%val,1)
      !stride = size(vec%val,2)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_extract_view


  subroutine phist_Dmvec_my_length(mvec_ptr, mylen, ierr) bind(C,name='phist_Dmvec_my_length_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    integer(C_INT32_T), intent(out) :: mylen
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      mylen = size(mvec%val,2)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_my_length


  subroutine phist_Dmvec_get_map(mvec_ptr, map_ptr, ierr) bind(C,name='phist_Dmvec_get_map_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    type(C_PTR),        intent(out) :: map_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      map_ptr = c_loc(mvec%map)
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_get_map


  subroutine phist_Dmvec_num_vectors(mvec_ptr, nvec, ierr) bind(C,name='phist_Dmvec_num_vectors_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value       :: mvec_ptr
    integer(C_INT),     intent(out) :: nvec
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      nvec = mvec%jmax-mvec%jmin+1
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_num_vectors


  subroutine phist_Dmvec_view_block(mvec_ptr, view_ptr, jmin, jmax, ierr) bind(C,name='phist_Dmvec_view_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    type(C_PTR),        intent(inout) :: view_ptr
    integer(C_INT),     value         :: jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, view
    !--------------------------------------------------------------------------------

    if( c_associated(mvec_ptr) ) then
      call c_f_pointer(mvec_ptr, mvec)
      if( c_associated(view_ptr) ) then
        call c_f_pointer(view_ptr,view)
        if( .not. view%is_view ) then
          write(*,*) 'mvec not a view!'
          flush(6)
          ierr = -88
          return
        end if
      else
        allocate(view)
        view_ptr = c_loc(view)
      end if
      view%is_view = .true.
      view%jmin = mvec%jmin+jmin
      view%jmax = mvec%jmin+jmax
      view%map=mvec%map
      view%val=>mvec%val
      ierr = 0
    else
      ierr = -88
    end if

  end subroutine phist_Dmvec_view_block


  subroutine phist_Dmvec_get_block(mvec_ptr, block_ptr, jmin, jmax, ierr) bind(C,name='phist_Dmvec_get_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr, block_ptr
    integer(C_INT),     value         :: jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, block
    type(MVec_t) :: view
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) .or. .not. c_associated(block_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)
    call c_f_pointer(block_ptr,block)

    ! create view and let mvec_add_mvec handle the rest!
    view%is_view = .true.
    view%map = mvec%map
    view%jmin = mvec%jmin+jmin
    view%jmax = mvec%jmin+jmax
    view%val=>mvec%val
    call mvec_add_mvec(1._8,view,0._8,block)

    ierr = 0

  end subroutine phist_Dmvec_get_block


  subroutine phist_Dmvec_set_block(mvec_ptr, block_ptr, jmin, jmax, ierr) bind(C,name='phist_Dmvec_set_block_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr, block_ptr
    integer(C_INT),     value         :: jmin, jmax
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, block
    type(MVec_t) :: view
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) .or. .not. c_associated(block_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)
    call c_f_pointer(block_ptr,block)

    ! create view and let mvec_add_mvec handle the rest!
    view%is_view = .true.
    view%map = mvec%map
    view%jmin = mvec%jmin+jmin
    view%jmax = mvec%jmin+jmax
    view%val=>mvec%val
    call mvec_add_mvec(1._8,block,0._8,view)

    ierr = 0

  end subroutine phist_Dmvec_set_block


  subroutine phist_Dmvec_gather_mvecs(mvec_ptr, block_ptr_list, nblocks, ierr) bind(C,name='phist_Dmvec_gather_mvecs_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    type(C_PTR),        intent(in)    :: block_ptr_list(*)
    integer(C_INT),     value         :: nblocks
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, tmp
    type(MVec_t), pointer :: block_list(:)
    integer :: i
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(mvec_ptr, mvec)

    allocate(block_list(nblocks))
    do i = 1, nblocks
      if( .not. c_associated(block_ptr_list(i)) ) then
        ierr = -88
        return
      end if
      call c_f_pointer(block_ptr_list(i),tmp)
      block_list(i)%jmin = tmp%jmin
      block_list(i)%jmax = tmp%jmax
      block_list(i)%is_view = .true.
      block_list(i)%map = tmp%map
      block_list(i)%val => tmp%val
    end do

    call mvec_gather_mvecs(mvec,block_list)
    deallocate(block_list)

    ierr = 0

  end subroutine phist_Dmvec_gather_mvecs


  subroutine phist_Dmvec_scatter_mvecs(mvec_ptr, block_ptr_list, nblocks, ierr) bind(C,name='phist_Dmvec_scatter_mvecs_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    type(C_PTR),        intent(in)    :: block_ptr_list(*)
    integer(C_INT),     value         :: nblocks
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec, tmp
    type(MVec_t), pointer :: block_list(:)
    integer :: i
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(mvec_ptr, mvec)

    allocate(block_list(nblocks))
    do i = 1, nblocks
      if( .not. c_associated(block_ptr_list(i)) ) then
        ierr = -88
        return
      end if
      call c_f_pointer(block_ptr_list(i),tmp)
      block_list(i)%jmin = tmp%jmin
      block_list(i)%jmax = tmp%jmax
      block_list(i)%is_view = .true.
      block_list(i)%map = tmp%map
      block_list(i)%val => tmp%val
    end do

    call mvec_scatter_mvecs(mvec,block_list)
    deallocate(block_list)

    ierr = 0

  end subroutine phist_Dmvec_scatter_mvecs


  subroutine phist_Dmvec_put_value(mvec_ptr, val, ierr) bind(C,name='phist_Dmvec_put_value_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    real(C_DOUBLE),     value         :: val
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    integer :: i
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

!$omp parallel do schedule(static)
    do i = 1, size(mvec%val,2), 1
      mvec%val(mvec%jmin:mvec%jmax,i) = val
    end do

    ierr = 0

  end subroutine phist_Dmvec_put_value


  subroutine phist_Dmvec_random(mvec_ptr, ierr) bind(C,name='phist_Dmvec_random_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

    call random_number(mvec%val(mvec%jmin:mvec%jmax,:))

    ierr = 0

  end subroutine phist_Dmvec_random


  subroutine phist_Dmvec_print(mvec_ptr, ierr) bind(C,name='phist_Dmvec_print_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    integer :: i
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

    do i = 1, size(mvec%val,2)
      write(*,*) mvec%val(mvec%jmin:mvec%jmax,i)
    end do
    flush(6)
    ierr = 0

  end subroutine phist_Dmvec_print


  subroutine phist_Dmvec_norm2(mvec_ptr, vnrm, ierr) bind(C,name='phist_Dmvec_norm2_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: mvec_ptr
    real(C_DOUBLE),     intent(out)   :: vnrm(*)
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: mvec
    integer :: nvec
    !--------------------------------------------------------------------------------

    if( .not. c_associated(mvec_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(mvec_ptr, mvec)

    nvec = mvec%jmax-mvec%jmin+1
    call mvec_norm2(mvec, vnrm(1:nvec))
    
    ierr = 0

  end subroutine phist_Dmvec_norm2


  subroutine phist_Dmvec_scale(x_ptr, alpha, ierr) bind(C,name='phist_Dmvec_scale_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: x_ptr
    real(C_DOUBLE),     value         :: alpha
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x
    !--------------------------------------------------------------------------------

    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)

    call mvec_scale(x,alpha)

    ierr = 0

  end subroutine phist_Dmvec_scale


  subroutine phist_Dmvec_vscale(x_ptr, alpha, ierr) bind(C,name='phist_Dmvec_vscale_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: x_ptr
    real(C_DOUBLE),     intent(in)    :: alpha(*)
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x
    !--------------------------------------------------------------------------------

    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)

    call mvec_vscale(x,alpha)

    ierr = 0

  end subroutine phist_Dmvec_vscale



  subroutine phist_Dmvec_add_mvec(alpha, x_ptr, beta, y_ptr, ierr) bind(C,name='phist_Dmvec_add_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    real(C_DOUBLE),     value         :: alpha, beta
    type(C_PTR),        value         :: x_ptr, y_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x, y
    !--------------------------------------------------------------------------------

    if( .not. c_associated(y_ptr) ) then
      ierr = -88
      return
    end if
    if( .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if

    call c_f_pointer(x_ptr, x)
    call c_f_pointer(y_ptr, y)

    call mvec_add_mvec(alpha, x, beta, y)
    ierr = 0

  end subroutine phist_Dmvec_add_mvec


  subroutine phist_Dmvec_vadd_mvec(alpha, x_ptr, beta, y_ptr, ierr) bind(C,name='phist_Dmvec_vadd_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    real(C_DOUBLE),     intent(in)    :: alpha(*)
    type(C_PTR),        value         :: x_ptr, y_ptr
    real(C_DOUBLE),     value         :: beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x, y
    integer :: nvec
    !--------------------------------------------------------------------------------

    if( beta .ne. 0 .and. .not. c_associated(y_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(y_ptr, y)
    nvec = y%jmax-y%jmin+1

    if( any(alpha(1:nvec) .ne. 0) .and. .not. c_associated(x_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(x_ptr, x)


    call mvec_vadd_mvec(alpha, x, beta, y)
    ierr = 0

  end subroutine phist_Dmvec_vadd_mvec


  subroutine phist_Dmvec_dot_mvec(x_ptr, y_ptr, dot, ierr) bind(C,name='phist_Dmvec_dot_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: x_ptr, y_ptr
    real(C_DOUBLE),     intent(out)   :: dot(*)
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: x, y
    !--------------------------------------------------------------------------------

    if( .not. c_associated(x_ptr) .or. .not. c_associated(y_ptr) ) then
      ierr = -88
      return
    end if
    call c_f_pointer(x_ptr, x)
    call c_f_pointer(y_ptr, y)

    call mvec_dot_mvec(x,y,dot)

    ierr = 0

  end subroutine phist_Dmvec_dot_mvec


  subroutine phist_Dmvec_times_sdMat(alpha, v_ptr, M_ptr, beta, w_ptr, ierr) bind(C,name='phist_Dmvec_times_sdMat_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    type(C_PTR),        value         :: M_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)

    call mvec_times_sdmat(alpha,v,M,beta,w)

    ierr = 0

  end subroutine phist_Dmvec_times_sdMat


  subroutine phist_Dmvec_times_sdMat_inplace(v_ptr, M_ptr, ierr) bind(C,name='phist_Dmvec_times_sdMat_inplace_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr
    type(C_PTR),        value         :: M_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v
    type(SDMat_t), pointer :: M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(M_ptr,M)

    call mvec_times_sdmat_inplace(v,M)

    ierr = 0

  end subroutine phist_Dmvec_times_sdMat_inplace


  subroutine phist_DmvecT_times_mvec(alpha, v_ptr, w_ptr, beta, M_ptr, ierr) bind(C,name='phist_DmvecT_times_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, w_ptr, M_ptr
    real(C_DOUBLE),     value         :: alpha, beta
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v, w
    type(SDMat_t), pointer :: M
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(w_ptr) .or. &
      & .not. c_associated(M_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(w_ptr,w)
    call c_f_pointer(M_ptr,M)

    call mvecT_times_mvec(alpha,v,w,beta,M)

    ierr = 0

  end subroutine phist_DmvecT_times_mvec


  subroutine phist_Dmvec_QR(v_ptr, R_ptr, ierr) bind(C,name='phist_Dmvec_QR_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),        value         :: v_ptr, R_ptr
    integer(C_INT),     intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(MVec_t), pointer :: v
    type(SDMat_t), pointer :: R
    !--------------------------------------------------------------------------------

    if( .not. c_associated(v_ptr) .or. &
      & .not. c_associated(R_ptr)      ) then
      ierr = -88
      return
    end if

    call c_f_pointer(v_ptr,v)
    call c_f_pointer(R_ptr,R)

    call mvec_QR(v,R,ierr)

  end subroutine phist_Dmvec_QR

end module mvec_module
