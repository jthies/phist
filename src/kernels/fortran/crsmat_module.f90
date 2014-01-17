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
  !> communication buffer for matrix vector multiplication of a matrix stored in crs format
  type CrsCommBuff_t
    integer :: nSendProcs                             !< local number of procs that need row blocks from this proc
    integer :: nRecvProcs                             !< local number of procs that store row blocks needed by this proc
    integer,      allocatable :: sendProcId(:)        !< mpi ids of procs data is sent to, has size (nSendProcs)
    integer,      allocatable :: recvProcId(:)        !< mpi ids of procs data is received from, has size (nRecvProcs)
    integer,      allocatable :: sendInd(:)           !< index in the send buffer of each proc, has size (nSendProcs+1)
    integer,      allocatable :: recvInd(:)           !< index in the receive buffer of each proc, has size (nRecvProcs+1)
    real(kind=8), allocatable :: sendData(:,:)        !< send buffer, has size (nb,sendBuffInd(nSendProcs+1)-1)
    real(kind=8), allocatable :: recvData(:,:)        !< recv buffer, has size (nb,recvBuffInd(nRecvProcs+1)-1)
    integer,      allocatable :: sendRowBlkInd(:)     !< local row block index of the data in the sendBuffer, has size (sendBuffInd(nSendProcs+1)-1)
    integer,      allocatable :: recvRowBlkInd(:)     !< global row block index of the data in the recvBuffer, has size (recvBuffInd(nRecvProcs+1)-1)
    integer,      allocatable :: sendRequests(:)      !< buffer for mpi send requests, has size (nSendProcs)
    integer,      allocatable :: recvRequests(:)      !< buffer for mpi receive requests, has size (nRecvProcs)
    integer,      allocatable :: sendStatus(:,:)      !< buffer for mpi send status, has size (MPI_STATUS_SIZE,nSendProcs)
    integer,      allocatable :: recvStatus(:,:)      !< buffer for mpi receive status, has size (MPI_STATUS_SIZE,nRecvProcs)
    integer,      allocatable :: recvIndices(:)       !< buffer for indices of received data, used for mpi_wait_some, has size(nRecvProcs)
  end type CrsCommBuff_t


  !==================================================================================
  !> distributed crs matrix
  type CrsMat_t
    !--------------------------------------------------------------------------------
    integer                      :: nRows                !< number of rows
    integer                      :: nCols                !< number of columns
    integer(kind=8)              :: nEntries             !< number of non-zero entries
    integer(kind=8), allocatable :: row_offset(:)        !< indices of rows in col_idx and val
    integer(kind=8), allocatable :: nonlocal_offset(:)   !< points to first nonlocal element in each row
    integer,         allocatable :: col_idx(:)           !< column indices
    real(kind=8),    allocatable :: val(:)               !< matrix entries
    type(Map_t)                  :: row_map
    type(CrsCommBuff_t)          :: comm_buff
    !--------------------------------------------------------------------------------
  end type CrsMat_t


contains


  !==================================================================================
  !> calculate the communication scheme that is needed for matrix vector
  !! multiplication for a distributed matrix in crs format.
  !! Also allocates necessary buffer space.
  subroutine setup_commBuff(mat,combuff)
    use mpi
    type(CrsMat_t),      intent(in)  :: mat
    type(CrsCommBuff_t), intent(out) :: combuff
    !--------------------------------------------------------------------------------
    integer,allocatable :: recvnum(:)
    integer,allocatable :: sendnum(:)
    integer(kind=8) :: i,j
    integer :: k, jProcIndex, jProc
    integer :: firstrow, lastrow
    integer :: buff_size, ierr
    logical,allocatable :: row_offset_counted(:)
    !--------------------------------------------------------------------------------

    ! determine nSendProcs and nRecvProcs
    allocate(recvnum(0:mat%row_map%nProcs-1))
    recvnum=0
    combuff%nRecvProcs=0
    firstrow = int(mat%row_map%distrib(mat%row_map%me))
    lastrow  = int(mat%row_map%distrib(mat%row_map%me+1)-1)
    do i=1,mat%nRows,1
      jProc=0
      do j=mat%row_offset(i),mat%row_offset(i+1)-1, 1
        if( mat%col_idx(j) .lt. firstrow .or. &
          & mat%col_idx(j) .gt. lastrow ) then
          do while( mat%col_idx(j) .ge. mat%row_map%distrib(jProc+1) )
            jProc=jProc+1
          end do
#ifdef TESTING
if( mat%col_idx(j) .lt. mat%row_map%distrib(jProc) ) then
  write(*,*) 'CRS sorting error! me', mat%row_map%me, 'idx', mat%col_idx(j), &
    &        'jProc', jProc, 'distrib', mat%row_map%distrib
  call exit(1)
end if
#endif
          if( recvnum(jProc) .eq. 0 ) then
            recvnum(jProc) = 1
            combuff%nRecvProcs=combuff%nRecvProcs+1
          end if
        end if
      end do
    end do
    allocate(sendnum(0:mat%row_map%nProcs-1))
    call mpi_allreduce(recvnum,sendnum,mat%row_map%nProcs,MPI_INTEGER,&
      &                MPI_SUM,mat%row_map%Comm,ierr)
    combuff%nSendProcs=sendnum(mat%row_map%me)
    deallocate(sendnum)

    ! allocate memory
    allocate(combuff%sendProcId(combuff%nSendProcs))
    allocate(combuff%recvProcId(combuff%nRecvProcs))
    allocate(combuff%sendInd(combuff%nSendProcs+1))
    allocate(combuff%recvInd(combuff%nRecvProcs+1))
    allocate(combuff%sendRequests(combuff%nSendProcs))
    combuff%sendRequests=MPI_REQUEST_NULL
    allocate(combuff%recvRequests(combuff%nRecvProcs))
    combuff%recvRequests=MPI_REQUEST_NULL
    allocate(combuff%sendStatus(MPI_STATUS_SIZE,combuff%nSendProcs))
    allocate(combuff%recvStatus(MPI_STATUS_SIZE,combuff%nRecvProcs))
    allocate(combuff%recvIndices(combuff%nRecvProcs))


    ! determine number of elements that need to be received from each proc
    allocate(row_offset_counted(mat%row_map%distrib(mat%row_map%nProcs)-1))
    combuff%recvInd(:) = 0
    row_offset_counted=.false.
    do i=1,mat%nRows,1
      jProc=-1
      jProcIndex=0
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        if( mat%col_idx(j) .lt. firstrow .or. &
          & mat%col_idx(j) .gt. lastrow ) then
          if( .not. row_offset_counted(mat%col_idx(j)) ) then
            row_offset_counted(mat%col_idx(j)) = .true.
            do while( mat%col_idx(j) .ge. mat%row_map%distrib(jProc+1) )
              jProc=jProc+1
              if( recvnum(jProc) .eq. 1 ) jProcIndex=jProcIndex+1
            end do
            combuff%recvInd(jProcIndex+1) = combuff%recvInd(jProcIndex+1)+1
            combuff%recvProcId(jProcIndex) = jProc
          end if
        end if
      end do
    end do

    ! determine number of elements that need to be sent to each proc
    do i=1,combuff%nRecvProcs,1
      jProc=combuff%recvProcId(i)
      call mpi_isend(combuff%recvInd(i+1),1,MPI_INTEGER,&
        &            jProc,10,mat%row_map%Comm,&
        &            combuff%recvRequests(i),ierr)
    end do


    combuff%sendInd(:)=0
    do i=1,combuff%nSendProcs,1
      call mpi_recv(combuff%sendInd(i+1),1,MPI_INTEGER,&
        &           MPI_ANY_SOURCE,10,mat%row_map%Comm,&
        &           combuff%sendStatus(:,i),ierr)
      combuff%sendProcId(i) = combuff%sendStatus(MPI_SOURCE,i)
    end do


    ! calculate total buffer size and allocate memory
    combuff%sendInd(1)=1
    do i=1,combuff%nSendProcs,1
      combuff%sendInd(i+1) = combuff%sendInd(i)+combuff%sendInd(i+1)
    end do
    buff_size=combuff%sendInd(combuff%nSendProcs+1)-1
    allocate(combuff%sendRowBlkInd(buff_size))
    call mpi_waitall(combuff%nRecvProcs,combuff%recvRequests,combuff%recvStatus,ierr)


    combuff%recvInd(1)=1
    do i=1,combuff%nRecvProcs,1
      combuff%recvInd(i+1) = combuff%recvInd(i)+combuff%recvInd(i+1)
    end do
    buff_size=combuff%recvInd(combuff%nRecvProcs+1)-1
    allocate(combuff%recvRowBlkInd(buff_size))


    ! determine row block indices of elements that we need to send/receive
    do i=1,combuff%nSendProcs,1
      j=combuff%sendInd(i)
      k=combuff%sendInd(i+1)
      call mpi_irecv(combuff%sendRowBlkInd(j:k-1),k-j,MPI_INTEGER,&
        &            combuff%sendProcId(i),20,mat%row_map%Comm,&
        &            combuff%sendRequests(i),ierr)
    end do

    row_offset_counted=.false.
    do i=1,mat%nRows,1
      jProc=-1
      jProcIndex=0
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        if( mat%col_idx(j) .lt. firstrow .or. &
          & mat%col_idx(j) .gt. lastrow ) then
          if( .not. row_offset_counted(mat%col_idx(j)) ) then
            row_offset_counted(mat%col_idx(j)) = .true.
            do while( mat%col_idx(j) .ge. mat%row_map%distrib(jProc+1) )
              jProc=jProc+1
              if( recvnum(jProc) .eq. 1 ) jProcIndex=jProcIndex+1
            end do
            combuff%recvRowBlkInd(combuff%recvInd(jProcIndex)) = mat%col_idx(j)
            combuff%recvInd(jProcIndex) = combuff%recvInd(jProcIndex)+1
          end if
        end if
      end do
    end do
    deallocate(row_offset_counted)
    deallocate(recvnum)
    
    ! restore recvInd
    do i=1,combuff%nRecvProcs-1,1
      combuff%recvInd(combuff%nRecvProcs-i+1)=combuff%recvInd(combuff%nRecvProcs-i)
    end do
    combuff%recvInd(1)=1
    
    ! sort and send data
    do i=1,combuff%nRecvProcs,1
      j=combuff%recvInd(i)
      k=combuff%recvInd(i+1)
      call sort(combuff%recvRowBlkInd(j:k-1))
      call mpi_isend(combuff%recvRowBlkInd(j:k-1),k-j,MPI_INTEGER,&
        &            combuff%recvProcId(i),20,mat%row_map%Comm,&
        &            combuff%recvRequests(i),ierr)
    end do

#ifdef TESTING
! check that recvRowBlkInd is globally sorted
do i = 2, size(combuff%recvRowBlkInd), 1
  if( combuff%recvRowBlkInd(i-1) .ge. combuff%recvRowBlkInd(i) ) then
    write(*,*) 'Error: recvRowBlkInd not sorted! me', mat%row_map%me, &
      &        'i', i, 'recvRowBlkInd', combuff%recvRowBlkInd
    flush(6)
    call exit(11)
  end if
end do
#endif

    call mpi_waitall(combuff%nSendProcs,combuff%sendRequests,&
      &              combuff%sendStatus,ierr)

   combuff%sendRowBlkInd = combuff%sendRowBlkInd - firstrow+1

    call mpi_waitall(combuff%nRecvProcs,combuff%recvRequests,&
      &              combuff%recvStatus,ierr)

  end subroutine setup_commBuff

  subroutine sort(array)
    use m_mrgrnk
    integer, intent(inout) :: array(1:)
    integer, allocatable :: tmp(:,:)
    integer :: i

    allocate(tmp(size(array),2))
    tmp(:,1) = array
    call mrgrnk(tmp(:,1),tmp(:,2))
    do i = 1, size(array)
      array(i) = tmp(tmp(i,2),1)
    end do

  end subroutine sort


  !> rearrange the elements of a matrix stored in crs format in a
  !! way that the elements are grouped by remote process (e.g. when performing
  !! a matrix vector multiplication, the elements that are multiplied with a
  !! subvector from a remote process are stored consecutivly in memory)
  subroutine sort_rows(mat)
    use m_mrgrnk
    type(CrsMat_t), intent(inout) :: mat

    integer :: i, n
    integer(kind=8) :: j, off
    integer :: firstrow, lastrow, recvBuffIndex
    integer, allocatable :: idx(:,:)
    real(kind=8), allocatable :: val(:)


    firstrow = mat%row_map%distrib(mat%row_map%me)
    lastrow  = mat%row_map%distrib(mat%row_map%me+1)-1

    ! distinguish between local elements and buffer indices
    allocate(mat%nonlocal_offset(mat%nRows))
    mat%nonlocal_offset(:) = mat%row_offset(1:mat%nRows)
    do i = 1, mat%nRows, 1
      do j = mat%row_offset(i), mat%row_offset(i+1)-1, 1
        if( mat%col_idx(j) .ge. firstrow .and. &
          & mat%col_idx(j) .le. lastrow        ) then
          ! local element, substract offset of this process
          mat%col_idx(j) = (mat%col_idx(j)-firstrow+1)
          ! increase nonlocal offset
          mat%nonlocal_offset(i) = mat%nonlocal_offset(i) + 1
        else
          ! nonlocal element, search its position in the comm buffer
          ! add nRows to sort them behind local elements
          mat%col_idx(j) = search_recvBuffIndex(mat%col_idx(j)) + mat%nRows
        end if
      end do
    end do

    allocate(idx(mat%nEntries,2))
    allocate(val(mat%nEntries))
    ! sort all rows
    do i = 1, mat%nRows
      off = mat%row_offset(i)
      n = mat%row_offset(i+1) - off
      if( n .gt. 0 ) then
        ! reuse idx for sorting
        idx(1:n,1) = mat%col_idx( off:off+n-1 )
        val(1:n)   = mat%val( off:off+n-1 )

        ! determine order of columns
        call mrgrnk(idx(1:n,1), idx(1:n,2))

        ! copy back sorted arrays
        do j = 1, n
          ! inverse sign, so local columns are positive, followed by negative indices in the recv buffer
          mat%col_idx( off+j-1 ) = idx(idx(j,2),1)
          mat%val( off+j-1 ) = val(idx(j,2))
        end do

      end if
      ! subtract mat%nRows from buffer indices, was only added for sorting!
      do j = mat%nonlocal_offset(i), mat%row_offset(i+1)-1, 1
        mat%col_idx(j) = mat%col_idx(j) - mat%nRows
      end do
    end do

  contains

  function search_recvBuffIndex(col) result(idx)
    integer, intent(in) :: col
    integer :: idx
    integer :: a, b, n, n_max


    a = 1
    b = size(mat%comm_buff%recvRowBlkInd)
    n_max = int(log(real(b))/log(2.))+1
    do n = 1, n_max, 1
      idx = (a+b)/2
      if( col .lt. mat%comm_buff%recvRowBlkInd(idx) ) then
        b = idx-1
      else if( col .gt. mat%comm_buff%recvRowBlkInd(idx) ) then
        a = idx+1
      else
        exit
      end if
    end do

    if( col .ne. mat%comm_buff%recvRowBlkInd(idx) ) then
      write(*,*) 'Error in search_recvBuffIndex! col:', col, 'a', a, 'b', b, 'idx', idx
      call exit(23)
    end if
  end function search_recvBuffIndex

  end subroutine sort_rows




  !==================================================================================
  !> multiply crsmat with mvec
  subroutine crsmat_times_mvec(alpha, A, shifts, x, beta, y)
    use, intrinsic :: iso_c_binding, only: C_INT
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),   intent(in)    :: alpha
    type(CrsMat_t), intent(inout) :: A
    real(kind=8),   intent(in)    :: shifts(*)
    type(MVec_t),   intent(in)    :: x
    real(kind=8),   intent(in)    :: beta
    type(MVec_t),   intent(inout) :: y
    !--------------------------------------------------------------------------------
    integer :: nvec, ldx, ldy, recvBuffSize
    logical :: strided_x, strided_y, strided
    logical :: y_is_aligned16, handled
    integer :: i, j, k, l, ierr
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
    strided = strided_x .or. strided_y
    nvec = x%jmax-x%jmin+1
    ldx = size(x%val,1)
    ldy = size(y%val,1)

    if( mod(loc(y%val(y%jmin,1)),16) .eq. 0 .and. (mod(ldy,2) .eq. 0 .or. ldy .eq. 1) ) then
      y_is_aligned16 = .true.
    else
      y_is_aligned16 = .false.
    end if

    write(*,*) 'spMVM with nvec =',nvec,', ldx =',ldx,', ldy =',ldy,', y_mem_aligned16 =', y_is_aligned16, 'on proc', A%row_map%me
    flush(6)


    ! exchange necessary elements
!write(*,*) 'CRS', A%row_map%me, 'sendRowBlkInd', A%comm_buff%sendRowBlkInd
!write(*,*) 'CRS', A%row_map%me, 'recvRowBlkInd', A%comm_buff%recvRowBlkInd

    ! start buffer irecv
    ! we could also set up persistent communication channels here... and use MPI_Startall later
    if( allocated(A%comm_buff%recvData) ) then
      if( size(A%comm_buff%recvData,1) .ne. nvec ) then
        deallocate(A%comm_buff%recvData)
        deallocate(A%comm_buff%sendData)
      end if
    end if
    if( .not. allocated(A%comm_buff%recvData) ) then
      allocate(A%comm_buff%recvData(nvec,A%comm_buff%recvInd(A%comm_buff%nRecvProcs+1)-1))
      allocate(A%comm_buff%sendData(nvec,A%comm_buff%sendInd(A%comm_buff%nSendProcs+1)-1))
    end if

    do i=1,A%comm_buff%nRecvProcs, 1
      k = A%comm_buff%recvInd(i)
      l = A%comm_buff%recvInd(i+1)
      call mpi_irecv(A%comm_buff%recvData(:,k:l-1),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
        &            A%comm_buff%recvProcId(i),3,A%row_map%Comm,A%comm_buff%recvRequests(i),ierr)
    end do

    ! start buffer isend
!$omp parallel
    do i=1,A%comm_buff%nSendProcs, 1
      k = A%comm_buff%sendInd(i)
      l = A%comm_buff%sendInd(i+1)
!$omp do
      do j = k, l-1, 1
        A%comm_buff%sendData(:,j) = x%val(x%jmin:x%jmax,A%comm_buff%sendRowBlkInd(j))
      end do
    end do
!$omp end parallel

    do i=1,A%comm_buff%nSendProcs, 1
      k = A%comm_buff%sendInd(i)
      l = A%comm_buff%sendInd(i+1)
      call mpi_isend(A%comm_buff%sendData(:,k:l-1),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
        &            A%comm_buff%sendProcId(i),3,A%row_map%Comm,A%comm_buff%sendRequests(i),ierr)
    end do

    if( A%comm_buff%nRecvProcs .gt. 0 ) then
      ! wait till all data arrived
      call mpi_waitall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,A%comm_buff%recvStatus,ierr)
    end if

    recvBuffSize = A%comm_buff%recvInd(A%comm_buff%nRecvProcs+1)-1


    handled = .false.
    !try to use NT stores if possible
    if( beta .eq. 0 .and. y_is_aligned16 ) then
      if( nvec .eq. 1 ) then
        if( .not. strided ) then
          call dspmvm_NT_1(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &              A%row_offset, A%nonlocal_offset, A%col_idx, A%val, shifts, x%val, &
            &              A%comm_buff%recvData, y%val)
          handled = .true.
        end if
      else if( nvec .eq. 2 ) then
        if( strided_x ) then
          call dspmvm_NT_strided_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &                      A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                      shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        else
          call dspmvm_NT_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &              A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &              shifts, x%val, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        end if
        handled = .true.
      else if( nvec .eq. 4 ) then
        if( strided_x ) then
          call dspmvm_NT_strided_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &                      A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                      shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        else
          call dspmvm_NT_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &              A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &              shifts, x%val, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        end if
        handled = .true.
      else if( nvec .eq. 8 ) then
        if( strided_x ) then
          call dspmvm_NT_strided_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &              A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                      shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        else
          call dspmvm_NT_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &              A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &              shifts, x%val, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        end if
        handled = .true.
      end if
    end if

    if( .not. handled ) then
      if( nvec .eq. 1 ) then
        if( strided ) then
          call dspmvm_strided_1(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &                   A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                   shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
        else
          call dspmvm_1(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &           A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &           shifts, x%val, A%comm_buff%recvData, beta, y%val)
        end if
        handled = .true.
      else if( nvec .eq. 2 ) then
        if( strided ) then
          call dspmvm_strided_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &                   A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                   shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
        else
          call dspmvm_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &           A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &           shifts, x%val, A%comm_buff%recvData, beta, y%val)
        end if
        handled = .true.
      else if( nvec .eq. 4 ) then
        if( strided ) then
          call dspmvm_strided_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &                   A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                   shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
        else
          call dspmvm_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &           A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &           shifts, x%val, A%comm_buff%recvData, beta, y%val)
        end if
        handled = .true.
      else if( nvec .eq. 8 ) then
        if( strided ) then
          call dspmvm_strided_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &                   A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &                   shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
        else
          call dspmvm_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
            &           A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
            &           shifts, x%val, A%comm_buff%recvData, beta, y%val)
        end if
        handled = .true.
      end if
    end if

!do i = 1, A%nRows, 1
  !write(*,*) 'CRS', A%row_map%me, 'localrow', i, 'entries', A%col_idx(A%row_offset(i):A%row_offset(i+1)-1)
!end do
!write(*,*) 'CRS', A%row_map%me, 'row_offset', A%row_offset
!write(*,*) 'CRS', A%row_map%me, 'halo_offset', A%nonlocal_offset
!write(*,*) 'CRS', A%row_map%me, 'col_idx', A%col_idx
!write(*,*) 'CRS', A%row_map%me, 'val', A%val
!write(*,*) 'CRS', A%row_map%me, 'alpha', alpha, 'beta', beta
!write(*,*) 'CRS', A%row_map%me, 'x', x%val(x%jmin:x%jmax,:)
!write(*,*) 'CRS', A%row_map%me, 'recvData', A%comm_buff%recvData
!write(*,*) 'CRS', A%row_map%me, 'y', y%val(y%jmin:y%jmax,:)
!flush(6)
    if( .not. handled ) then
      call dspmvm_generic(nvec, A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
        &                 A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
        &                 shifts, x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
    end if


    if( A%comm_buff%nSendProcs .gt. 0 ) then
      ! just make sure the buffers are not used any more...
      call mpi_waitall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
    end if

    !--------------------------------------------------------------------------------
  end subroutine crsmat_times_mvec


  !==================================================================================
  !> read MatrixMarket file
  subroutine phist_DcrsMat_read_mm(A_ptr, filename_len, filename_ptr, ierr) bind(C,name='phist_DcrsMat_read_mm_f')
    use, intrinsic :: iso_c_binding
    use env_module, only: newunit
    use mpi
    use m_mrgrnk
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
    real(kind=8), allocatable :: val(:)
    integer :: i, j, off, n, globalRows, globalCols, globalEntries
    !--------------------------------------------------------------------------------

    do i = 1, filename_len
      filename(i:i) = filename_ptr(i)
    end do

    ! open the file
    write(*,*) 'reading file:', filename
    flush(6)
    open(unit   = newunit(funit), file    = filename, &
      &  action = 'read',         status  = 'old',    &
      &  iostat = ierr)
    if( ierr .ne. 0 ) return

    ! read first line
    read(funit,'(A)') line
    write(*,*) line
    flush(6)
    if( trim(line) .ne. '%%MatrixMarket matrix coordinate real general' ) then
      write(*,*) 'unsupported format'
      flush(6)
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
    read(funit,*) globalRows, globalCols, globalEntries
    write(*,*) 'CrsMat:', globalRows, globalCols, globalEntries
    flush(6)

    call map_setup(A%row_map, MPI_COMM_WORLD, int(globalRows,kind=8), ierr)
    if( ierr .ne. 0 ) return

    A%nRows = A%row_map%nlocal(A%row_map%me)
    A%nCols = A%nRows

    ! allocate temporary buffers
    allocate(idx(globalEntries,2))
    allocate(val(globalEntries))

    ! read data
    j = 0
    do i = 1, globalEntries, 1
      read(funit,*) idx(j+1,1), idx(j+1,2), val(j+1)
      if( idx(j+1,1) .ge. A%row_map%distrib(A%row_map%me) .and. &
        & idx(j+1,1) .lt. A%row_map%distrib(A%row_map%me+1)     ) then
        j = j + 1
        ! already subtract offset of this proc
        idx(j,1) = idx(j,1) -  A%row_map%distrib(A%row_map%me) + 1
      end if
    end do
    A%nEntries = j

    ! close the file
    close(funit)


    ! allocate crs matrix
    allocate(A%row_offset(A%nRows+1))
    allocate(A%col_idx(A%nEntries))
    allocate(A%val(A%nEntries))

    ! try to respect NUMA
!$omp workshare
    A%row_offset = 0
    A%col_idx = 0
    A%val = 0._8
!$omp end workshare

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
          A%col_idx( off+j-1 ) = idx(idx(j,2),1)
          A%val( off+j-1 ) = val(idx(j,2))
        end do
      end if
    end do

    call setup_commBuff(A, A%comm_buff)
    call sort_rows(A)

    write(*,*) 'created new crsMat with dimensions', A%nRows, A%nCols, A%nEntries
    flush(6)
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

#ifdef TESTING
    write(*,*) 'deleting crsMat at address', A_ptr
    flush(6)
#endif
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
    real(kind=8), allocatable :: shifts(:)
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

    allocate(shifts(x%jmax-x%jmin+1))
    shifts = 0._8

    call crsmat_times_mvec(alpha,A,shifts,x,beta,y)

    ierr = 0
  end subroutine phist_DcrsMat_times_mvec


  !==================================================================================
  subroutine phist_DcrsMat_times_mvec_vadd_mvec(alpha, A_ptr, shifts, x_ptr, beta, y_ptr, ierr) bind(C,name='phist_DcrsMat_times_mvec_vadd_mvec_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    real(C_DOUBLE),   value         :: alpha, beta
    real(C_DOUBLE),   intent(in)    :: shifts(*)
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

    call crsmat_times_mvec(alpha,A,shifts,x,beta,y)

    ierr = 0
  end subroutine phist_DcrsMat_times_mvec_vadd_mvec


end module crsmat_module

