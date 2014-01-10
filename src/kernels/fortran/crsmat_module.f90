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


  !==================================================================================
  !> offsets and indices for proc-wise sorted bcrs
  !! only used internally!
  type ProcWiseIndices_t
    integer :: nd                                     !< number of rows with elements from this proc
    integer :: offset                                 !< row block offset of first row from this proc
    integer(kind=8), allocatable :: row_offset(:)     !< index of elements in this proc, has size (nd+1)
  end type ProcWiseIndices_t


contains

  !==================================================================================
  !> rearranges the entries in CrsMat_t%val and CrsMat_t%col_idx in a given order
  pure subroutine reordCRS(mat, new_ind)
    type(CrsMat_t),  intent(inout) :: mat
    integer(kind=8), intent(inout) :: new_ind(:)

    integer(kind=8) :: j,l
    integer :: buff_col_idx
    integer(kind=8) :: buff_new_ind
    real(kind=8) :: buff_value

    ! now sort the data in place (we have already the future position,
    ! so this can be done by O(n) swaps)
    do j=1,mat%nEntries, 1
      do while( new_ind(j) .ne. j )
        l = new_ind(j)
        ! swap j and l

        buff_value = mat%val(l)
        buff_col_idx = mat%col_idx(l)
        buff_new_ind = new_ind(l)

        mat%val(l) = mat%val(j)
        mat%col_idx(l) = mat%col_idx(j)
        new_ind(l) = new_ind(j)

        mat%val(j) = buff_value
        mat%col_idx(j) = buff_col_idx
        new_ind(j) = buff_new_ind

      end do
    end do

  end subroutine reordCRS


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
    
    call mpi_waitall(combuff%nSendProcs,combuff%sendRequests,&
      &              combuff%sendStatus,ierr)

   combuff%sendRowBlkInd = combuff%sendRowBlkInd - firstrow+1

    call mpi_waitall(combuff%nRecvProcs,combuff%recvRequests,&
      &              combuff%recvStatus,ierr)

  contains

    subroutine sort(array)
      use m_mrgrnk
      integer, intent(inout) :: array(1:)
      integer, allocatable :: tmp(:,:)

      allocate(tmp(size(array),2))
      tmp(:,1) = array
      call mrgrnk(tmp(:,1),tmp(:,2))
      do i = 1, size(array)
        array(tmp(i,2)) = tmp(i,1)
      end do
    end subroutine sort
  end subroutine setup_commBuff


  !> rearrange the elements of a matrix stored in crs format in a
  !! way that the elements are grouped by remote process (e.g. when performing
  !! a matrix vector multiplication, the elements that are multiplied with a
  !! subvector from a remote process are stored consecutivly in memory)
  subroutine sortlnlcrs(mat)
    type(CrsMat_t), intent(inout) :: mat

    integer :: i, k, jProcIndex
    integer(kind=8) :: j, tmp
    integer(kind=8),allocatable :: proc_ind(:)
    integer :: firstrow, lastrow, recvBuffIndex
    type(ProcWiseIndices_t), allocatable :: lnlInd(:)


    firstrow = mat%row_map%distrib(mat%row_map%me)
    lastrow  = mat%row_map%distrib(mat%row_map%me+1)-1

    ! first create new row_offset indices per proc
    allocate(lnlInd(0:mat%comm_buff%nRecvProcs))
    allocate(proc_ind(mat%nEntries))
    ! give each element a proc_ind and get row block offset
    ! and number of rows per proc
    do i=0,mat%comm_buff%nRecvProcs, 1
      lnlInd(i)%offset=mat%nRows
      lnlInd(i)%nd=0
    end do
    do i=1,mat%nRows, 1
      jProcIndex=1
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        if( mat%col_idx(j) .ge. firstrow .and. &
          & mat%col_idx(j) .le. lastrow ) then
          ! is a local element
          proc_ind(j) = 0
        else
          ! non local element
          do while(mat%col_idx(j) .ge. &
            &      mat%row_map%distrib(mat%comm_buff%recvProcId(jProcIndex)+1) )
            jProcIndex=jProcIndex+1
          end do
          proc_ind(j) = jProcIndex
        end if

        lnlInd(proc_ind(j))%offset = min(i,lnlInd(proc_ind(j))%offset)
        lnlInd(proc_ind(j))%nd = i-lnlInd(proc_ind(j))%offset+1
      end do
    end do

    ! allocate memory for indices
    do i=0,mat%comm_buff%nRecvProcs, 1
      allocate(lnlInd(i)%row_offset( lnlInd(i)%nd+1 ))
      lnlInd(i)%row_offset=0
    end do

    ! go through all elements and count the elements per proc
    do i=1,mat%nRows, 1
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        lnlInd(proc_ind(j))%row_offset(i-lnlInd(proc_ind(j))%offset+2)=& 
        &lnlInd(proc_ind(j))%row_offset(i-lnlInd(proc_ind(j))%offset+2)+1
      end do
    end do

    ! make indices from elem counts
    j=1
    do i=0,mat%comm_buff%nRecvProcs, 1
      do k=1,lnlInd(i)%nd+1, 1
        j=j+lnlInd(i)%row_offset(k)
        lnlInd(i)%row_offset(k)=j
      end do
    end do



    ! we have created all necessary indices, now we need to sort the
    ! data in the matrix
    ! save the new position in proc_ind
    do i=1,mat%nRows,1
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        tmp = lnlInd(proc_ind(j))%& 
          &           row_offset(i-lnlInd(proc_ind(j))%offset+1)
        lnlInd(proc_ind(j))%& 
          &           row_offset(i-lnlInd(proc_ind(j))%offset+1) = &
          &           tmp + 1
        proc_ind(j) = tmp
      end do
    end do

    ! restore indices
    do i=0,mat%comm_buff%nRecvProcs, 1
      do k=lnlInd(i)%nd,1,-1
        lnlInd(i)%row_offset(k+1)=lnlInd(i)%row_offset(k)
      end do
      if( i .eq. 0 ) then
        lnlInd(i)%row_offset(1)=1
      else
        lnlInd(i)%row_offset(1)=lnlInd(i-1)%row_offset(&
          &                 lnlInd(i-1)%nd+1)
      end if
    end do

    ! reorder the elements
    call reordCRS(mat, proc_ind)
    deallocate(proc_ind)

    ! make the mat%col_idx local
    do i=1,lnlInd(0)%nd,1
      do j=lnlInd(0)%row_offset(i),lnlInd(0)%row_offset(i+1)-1,1
        mat%col_idx(j)=mat%col_idx(j)-mat%row_map%distrib(mat%row_map%me)+1
      end do
    end do

    do k=1,mat%comm_buff%nRecvProcs, 1
      do i=1,lnlInd(k)%nd,1
        recvBuffIndex = mat%comm_buff%recvInd(k)
        do j=lnlInd(k)%row_offset(i),lnlInd(k)%row_offset(i+1)-1,1
          do while( mat%col_idx(j) .gt. mat%comm_buff%recvRowBlkInd(recvBuffIndex) )
            recvBuffIndex=recvBuffIndex+1
          end do
          mat%col_idx(j) = -recvBuffIndex
        end do
      end do
    end do

    ! get pointer to first nonlocal element in each row
    allocate(mat%nonlocal_offset(mat%nRows))
    do i=1,mat%nRows,1
      do j = mat%row_offset(i), mat%row_offset(i+1)-1, 1
        if( mat%col_idx(j) .lt. 0 ) exit
      end do
      mat%nonlocal_offset(i) = j
    end do

  end subroutine sortlnlcrs




  !==================================================================================
  !> multiply crsmat with mvec
  subroutine crsmat_times_mvec(alpha, A, x, beta, y)
    use, intrinsic :: iso_c_binding, only: C_INT
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),   intent(in)    :: alpha
    type(CrsMat_t), intent(inout) :: A
    type(MVec_t),   intent(in)    :: x
    real(kind=8),   intent(in)    :: beta
    type(MVec_t),   intent(inout) :: y
    !--------------------------------------------------------------------------------
    integer :: nvec, ldx, ldy, recvBuffSize
    logical :: strided_x, strided_y
    logical :: y_is_aligned16
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
    nvec = x%jmax-x%jmin+1
    ldx = size(x%val,1)
    ldy = size(y%val,1)

    if( mod(loc(y%val(y%jmin,1)),16) .eq. 0 .and. mod(ldy,2) .eq. 0 ) then
      y_is_aligned16 = .true.
    else
      y_is_aligned16 = .false.
    end if

    write(*,*) 'spMVM with nvec =',nvec,', ldx =',ldx,', ldy =',ldy,', y_mem_aligned16 =', y_is_aligned16
    flush(6)


    ! exchange necessary elements

    ! start buffer irecv
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
    do i=1,A%comm_buff%nSendProcs, 1
      k = A%comm_buff%sendInd(i)
      l = A%comm_buff%sendInd(i+1)
      do j = k, l-1, 1
        A%comm_buff%sendData(:,j) = x%val(x%jmin:x%jmax,A%comm_buff%sendRowBlkInd(j))
      end do
      call mpi_isend(A%comm_buff%sendData(:,k:l-1),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
        &            A%comm_buff%sendProcId(i),3,A%row_map%Comm,A%comm_buff%sendRequests(i),ierr)
    end do

    ! wait till all data arrived
    call mpi_waitall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,A%comm_buff%recvStatus,ierr)

    recvBuffSize = A%comm_buff%recvInd(A%comm_buff%nRecvProcs+1)-1


    ! try to use NT stores if possible
    !if( beta .eq. 0 .and. y_is_aligned16 ) then
      !if( nvec .eq. 1 ) then
        !if( .not. strided_x ) then
          !call dspmvm_NT_1(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, x%val, &
            !&              A%comm_buff%recvData, y%val)
          !return
        !end if
      !else if( nvec .eq. 2 ) then
        !if( strided_x ) then
          !call dspmvm_NT_strided_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            !&                      x%val(x%jmin,1), ldx, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        !else
          !call dspmvm_NT_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            !&              x%val, y%val(y%jmin,1), ldy)
        !end if
        !return
      !else if( nvec .eq. 4 ) then
        !if( strided_x ) then
          !call dspmvm_NT_strided_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            !&                      x%val(x%jmin,1), ldx, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        !else
          !call dspmvm_NT_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            !&              x%val, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        !end if
        !return
      !else if( nvec .eq. 8 ) then
        !if( strided_x ) then
          !call dspmvm_NT_strided_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            !&                      x%val(x%jmin,1), ldx, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        !else
          !call dspmvm_NT_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
            !&              x%val, A%comm_buff%recvData, y%val(y%jmin,1), ldy)
        !end if
        !return
      !end if
    !end if

    !if( nvec .eq. 1 ) then
      !if( strided_x .or. strided_y ) then
        !call dspmvm_strided_1(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&                   x%val(x%jmin,1), ldx, beta, y%val(y%jmin,1), ldy)
      !else
        !call dspmvm_1(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&           x%val, A%comm_buff%recvData, beta, y%val)
      !end if
      !return
    !else if( nvec .eq. 2 ) then
      !if( strided_x .or. strided_y ) then
        !call dspmvm_strided_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&                   x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
      !else
        !call dspmvm_2(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&           x%val, A%comm_buff%recvData, beta, y%val)
      !end if
      !return
    !else if( nvec .eq. 4 ) then
      !if( strided_x .or. strided_y ) then
        !call dspmvm_strided_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&                   x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
      !else
        !call dspmvm_4(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&           x%val, A%comm_buff%recvData, beta, y%val)
      !end if
      !return
    !else if( nvec .eq. 8 ) then
      !if( strided_x .or. strided_y ) then
        !call dspmvm_strided_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&                   x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)
      !else
        !call dspmvm_8(A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, A%row_offset, A%col_idx, A%val, &
          !&           x%val, A%comm_buff%recvData, beta, y%val)
      !end if
      !return
    !end if

    call dspmvm_generic(nvec, A%nrows, recvBuffSize, A%ncols, A%nEntries, alpha, &
      &                 A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
      &                 x%val(x%jmin,1), ldx, A%comm_buff%recvData, beta, y%val(y%jmin,1), ldy)


    ! just make sure the buffers are not used any more...
    call mpi_waitall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)

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
    read(funit,*) A%nRows, A%nCols, A%nEntries
    write(*,*) 'CrsMat:', A%nRows, A%nCols, A%nEntries
    flush(6)

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
          A%col_idx( off+idx(j,2)-1 ) = idx(j,1)
          A%val( off+idx(j,2)-1 ) = val(j)
        end do
      end if
    end do

    call setup_commBuff(A, A%comm_buff)
    call sortlnlCrs(A)

    if( A%row_map%me .eq. 0 ) then
      write(*,*) 'created new crsMat with dimensions', A%nRows, A%nCols, A%nEntries
      flush(6)
    end if
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

