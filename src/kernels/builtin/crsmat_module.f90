!> \file crsmat_module.f90
!! defines crsmat_module, the sparseMat implementation of phist builtin kernels
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>"
!! \author "Jonas Thies <Jonas.Thies@DLR.de>"
!!

#include "phist_config.h"
#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#endif

!> Implements phist_DsparseMat_* routines for the builtin kernels
module crsmat_module
  use, intrinsic :: iso_c_binding
  use map_module, only: Map_t, map_setup, map_compatible_map
  use mvec_module, only: MVec_t, mvec_scale
  implicit none

  private


#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
#define G_GIDX_T C_INT64_T
#else
#ifdef GHOST_HAVE_LONGIDX_LOCAL
#define G_LIDX_T C_INT64_T
#else
#define G_LIDX_T C_INT32_T
#endif
#ifdef GHOST_HAVE_LONGIDX_GLOBAL
#define G_GIDX_T C_INT64_T
#else
#define G_GIDX_T C_INT32_T
#endif
#endif

! this macro decides how many mvecs can be communicated at a time
! (should be >=1 for spMVM and >=2 for the CARP kernel with complex shifts)
#define NBLOCKS 2

  public :: CrsMat_t
  !public :: phist_DcrsMat_read_mm
  !public :: phist_DcrsMat_create_fromRowFunc
  !public :: phist_DcrsMat_delete
  !public :: phist_DcrsMat_times_mvec
  !public :: phist_DcrsMat_times_mvec_vadd_mvec
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
    ! buffers to allow sending several (nb') blocks of nb vectors each. nb' is set by the macro NBLOCKS above
    real(kind=8), allocatable :: sendData(:,:,:)        !< send buffer, has size (nb,sendBuffInd(nSendProcs+1)-1,nb')
    real(kind=8), allocatable :: recvData(:,:,:)        !< recv buffer, has size (nb,recvBuffInd(nRecvProcs+1)-1,nb')
    integer,      allocatable :: sendRowBlkInd(:)     !< local row block index of the data in the sendBuffer, has size (sendBuffInd(nSendProcs+1)-1)
    integer,      allocatable :: sendBuffInd(:,:)     !< sorted Variant of sendRowBlkInd, hast size (size(sendRowBlkInd),2), first entry is the local row block index, second is the index in the send buffer
    integer(kind=8), allocatable :: recvRowBlkInd(:)     !< global row block index of the data in the recvBuffer, has size (recvBuffInd(nRecvProcs+1)-1)
    integer,      allocatable :: sendRequests(:)      !< buffer for mpi send requests, has     size (nSendProcs,nb')
    integer,      allocatable :: recvRequests(:)      !< buffer for mpi receive requests, has size (nRecvProcs,nb')
    integer,      allocatable :: inv_sendRequests(:)  !< MPI requests for reversed communication in CARP-CG
    integer,      allocatable :: inv_recvRequests(:)  !< MPI requests for reversed communication in CARP-CG
    integer,      allocatable :: sendStatus(:,:)      !< buffer for mpi send status, has size (MPI_STATUS_SIZE,nSendProcs,nb')
    integer,      allocatable :: recvStatus(:,:)      !< buffer for mpi receive status, has size (MPI_STATUS_SIZE,nRecvProcs,nb')
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
    integer(kind=8), allocatable :: global_col_idx(:)    !< global column indices
    integer(kind=8), allocatable :: global_row_idx(:)    !< from possible reordering: original global indices of local rows
    type(Map_t)                  :: row_map
    type(CrsCommBuff_t)          :: comm_buff
    !--------------------------------------------------------------------------------
  end type CrsMat_t

  !! a simple sparse vector object, used to implement averaging 
  !! of received values in carp_sweep without tedium.           
  !! TODO: this object may cause NUMA problems
  type sparseArray_t
  
  integer :: n !! maximum idx allowed
  integer :: nnz !! number of elements stored
  integer, allocatable :: idx(:) !! my local indices into a full vector representation
  real(kind=8), allocatable :: val(:) !! val(i)=v(idx(i))
  
  end type sparseArray_t

  !> interface of function-ptr for crsMat_create_fromRowFunc
  abstract interface
    subroutine matRowFunc(row, nnz, cols, vals)
      use, intrinsic :: iso_c_binding
      integer(G_GIDX_T), value :: row
      integer(G_LIDX_T), intent(inout) :: nnz
      integer(G_GIDX_T), intent(inout) :: cols(*)
      real(C_DOUBLE),    intent(inout) :: vals(*)
    end subroutine matRowFunc
  end interface

contains


  !==================================================================================
  !> calculate the communication scheme that is needed for matrix vector
  !! multiplication for a distributed matrix in crs format.
  !! Also allocates necessary buffer space.
  subroutine setup_commBuff(mat,combuff)
    use mpi
    type(CrsMat_t),      intent(in)    :: mat
    type(CrsCommBuff_t), intent(inout) :: combuff
    !--------------------------------------------------------------------------------
    integer,allocatable :: recvnum(:)
    integer,allocatable :: sendnum(:)
    integer(kind=8) :: i,j
    integer :: k, jProcIndex, jProc
    integer(kind=8) :: firstrow, lastrow
    integer :: buff_size, ierr
    logical,allocatable :: row_offset_counted(:)
    integer(kind=8), allocatable :: sendRowBlkInd(:)
    !--------------------------------------------------------------------------------

    ! determine nSendProcs and nRecvProcs
    allocate(recvnum(0:mat%row_map%nProcs-1))
    recvnum=0
    combuff%nRecvProcs=0
    firstrow = mat%row_map%distrib(mat%row_map%me)
    lastrow  = mat%row_map%distrib(mat%row_map%me+1)-1
    do i=1,mat%nRows,1
      jProc=0
      do j=mat%row_offset(i),mat%row_offset(i+1)-1, 1
        if( mat%global_col_idx(j) .lt. firstrow .or. &
          & mat%global_col_idx(j) .gt. lastrow ) then
          do while( mat%global_col_idx(j) .ge. mat%row_map%distrib(jProc+1) )
            jProc=jProc+1
          end do
#if defined(TESTING) && PHIST_OUTLEV >= 4
if( mat%global_col_idx(j) .lt. mat%row_map%distrib(jProc) ) then
  write(*,*) 'CRS sorting error! me', mat%row_map%me, 'idx', mat%global_col_idx(j), &
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
    allocate(combuff%sendRequests(combuff%nSendProcs*NBLOCKS))
    combuff%sendRequests=MPI_REQUEST_NULL
    allocate(combuff%recvRequests(combuff%nRecvProcs*NBLOCKS))
    combuff%recvRequests=MPI_REQUEST_NULL
    allocate(combuff%inv_recvRequests(combuff%nSendProcs*NBLOCKS))
    combuff%inv_recvRequests=MPI_REQUEST_NULL
    allocate(combuff%inv_sendRequests(combuff%nRecvProcs*NBLOCKS))
    combuff%inv_sendRequests=MPI_REQUEST_NULL
    allocate(combuff%sendStatus(MPI_STATUS_SIZE,combuff%nSendProcs*NBLOCKS))
    allocate(combuff%recvStatus(MPI_STATUS_SIZE,combuff%nRecvProcs*NBLOCKS))
    allocate(combuff%recvIndices(combuff%nRecvProcs))


    ! determine number of elements that need to be received from each proc
    allocate(row_offset_counted(mat%row_map%distrib(mat%row_map%nProcs)-1))
    combuff%recvInd(:) = 0
    row_offset_counted=.false.
    do i=1,mat%nRows,1
      jProc=-1
      jProcIndex=0
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        if( mat%global_col_idx(j) .lt. firstrow .or. &
          & mat%global_col_idx(j) .gt. lastrow ) then
          if( .not. row_offset_counted(mat%global_col_idx(j)) ) then
            row_offset_counted(mat%global_col_idx(j)) = .true.
            do while( mat%global_col_idx(j) .ge. mat%row_map%distrib(jProc+1) )
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
    allocate(sendRowBlkInd(buff_size))
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
      call mpi_irecv(sendRowBlkInd(j:k-1),int(k-j),MPI_INTEGER8,&
        &            combuff%sendProcId(i),20,mat%row_map%Comm,&
        &            combuff%sendRequests(i),ierr)
    end do

    row_offset_counted=.false.
    do i=1,mat%nRows,1
      jProc=-1
      jProcIndex=0
      do j=mat%row_offset(i),mat%row_offset(i+1)-1,1
        if( mat%global_col_idx(j) .lt. firstrow .or. &
          & mat%global_col_idx(j) .gt. lastrow ) then
          if( .not. row_offset_counted(mat%global_col_idx(j)) ) then
            row_offset_counted(mat%global_col_idx(j)) = .true.
            do while( mat%global_col_idx(j) .ge. mat%row_map%distrib(jProc+1) )
              jProc=jProc+1
              if( recvnum(jProc) .eq. 1 ) jProcIndex=jProcIndex+1
            end do
            combuff%recvRowBlkInd(combuff%recvInd(jProcIndex)) = mat%global_col_idx(j)
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
      call mpi_isend(combuff%recvRowBlkInd(j:k-1),int(k-j),MPI_INTEGER8,&
        &            combuff%recvProcId(i),20,mat%row_map%Comm,&
        &            combuff%recvRequests(i),ierr)
    end do

#if defined(TESTING) && PHIST_OUTLEV >= 4
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


    combuff%sendRowBlkInd = int(sendRowBlkInd - firstrow+1,kind=4)

    call mpi_waitall(combuff%nRecvProcs,combuff%recvRequests,&
      &              combuff%recvStatus,ierr)


    ! get sorted variant of sendRowBlkInd for faster copy-to-buffer later
    allocate(combuff%sendBuffInd(size(combuff%sendRowBlkInd),2))
    do i = 1, size(combuff%sendRowBlkInd), 1
      combuff%sendBuffInd(i,1) = combuff%sendRowBlkInd(i)
      combuff%sendBuffInd(i,2) = int(i,kind=4)
    end do
    call sortInd(combuff%sendBuffInd(:,1), combuff%sendBuffInd(:,2))
    !write(*,*) 'sendBuffInd 1', combuff%sendBuffInd(:,1)
    !write(*,*) 'sendBuffInd 2', combuff%sendBuffInd(:,2)

  end subroutine setup_commBuff

  !! these private subroutines need to be caleed before
  !! an spMVM if nvec changes (not by the user, of course)
  subroutine delete_buffers(B)
  USE mpi
  implicit none
  TYPE(CrsCommBuff_t) :: B
  integer :: i, j, ierr

  ! deallocate old receive buffer
  do j = 1, NBLOCKS
  do i = 1, B%nRecvProcs
     call mpi_request_free(B%recvRequests((j-1)*B%nRecvProcs+i), ierr)
  end do
  do i = 1, B%nSendProcs
     call mpi_request_free(B%sendRequests((j-1)*B%nSendProcs+i), ierr)
  end do
  end do
  deallocate(B%recvData)
  deallocate(B%sendData)
  
  end subroutine delete_buffers

  !! allocate buffers for sending/receiving vector entries,
  !! and setup persistent communication. nvec is nb above,
  !! and nblock is nb' above.
  subroutine alloc_buffers(B,comm,nvec,ierr)
  USE mpi  
  implicit none
  
  TYPE(CrsCommBuff_t) :: B
  INTEGER, intent(in) :: comm
  INTEGER, intent(in) :: nvec
  INTEGER, intent(out) :: ierr
  INTEGER :: i,j,k,l
  INTEGER(KIND=8) :: recvBuffSize,sendBuffSize

  ! we allocate at least one element so that we can pass the first element of the
  ! buffer to some subroutine (or the first element of the trailing dimension subarray)

  ! allocate new receive buffer
  recvBuffSize=MAX(B%recvInd(B%nRecvProcs+1)-1,1)
  allocate(B%recvData(nvec,recvBuffSize,NBLOCKS),stat=ierr)
  if (ierr/=0) return

  ! allocate new send buffer
  sendBuffSize=MAX(B%sendInd(B%nSendProcs+1)-1,1)
  allocate(B%sendData(nvec,sendBuffSize,NBLOCKS),stat=ierr)
  if (ierr/=0) return

  ! setup persistent communication for spMVM
  do j = 1, NBLOCKS
  do i = 1, B%nRecvProcs
    k = B%recvInd(i)
    l = B%recvInd(i+1)
    call mpi_recv_init(B%recvData(1,k,j),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
      &        B%recvProcId(i),3,comm,B%recvRequests((j-1)*B%nRecvProcs+i),ierr)
    if (ierr/=0) return
  end do

  do i=1,B%nSendProcs, 1
    k = B%sendInd(i)
    l = B%sendInd(i+1)
    call mpi_ssend_init(B%sendData(1,k,j),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
      &                 B%sendProcId(i),3,comm,B%sendRequests((j-1)*B%nSendProcs+i),ierr)
    if (ierr/=0) return
  end do
  end do

  ! setup persistent communication inverse data transfers
  ! (sending from recvBuf and vice versa, used in CARP)
  do j = 1, NBLOCKS
  do i = 1, B%nRecvProcs
    k = B%recvInd(i)
    l = B%recvInd(i+1)
    call mpi_send_init(B%recvData(1,k,j),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
      &              B%recvProcId(i),4,comm,B%inv_sendRequests((j-1)*B%nRecvProcs+i),ierr)
    if (ierr/=0) return
  end do

  do i=1,B%nSendProcs, 1
    k = B%sendInd(i)
    l = B%sendInd(i+1)
    call mpi_recv_init(B%sendData(1,k,j),(l-k)*nvec,MPI_DOUBLE_PRECISION,&
      &                B%sendProcId(i),4,comm,B%inv_recvRequests((j-1)*B%nSendProcs+i),ierr)
    if (ierr/=0) return
  end do
  end do

  end subroutine alloc_buffers

  subroutine sort(array)
    use m_mrgrnk
    integer(kind=8), intent(inout) :: array(1:)
    integer(kind=8), allocatable :: tmp(:,:)
    integer :: i

    allocate(tmp(size(array),2))
    tmp(:,1) = array
    call mrgrnk(tmp(:,1),tmp(:,2))
    do i = 1, size(array)
      array(i) = tmp(tmp(i,2),1)
    end do

  end subroutine sort

  subroutine sortInd(array, ind)
    use m_mrgrnk
    integer(kind=4), intent(inout) :: array(1:)
    integer(kind=4), intent(out) :: ind(1:)
    integer(kind=8), allocatable :: tmp(:,:)
    integer :: i

    allocate(tmp(size(array),2))
    tmp(:,1) = int(array,kind=8)
    call mrgrnk(tmp(:,1),tmp(:,2))
    do i = 1, size(array)
      array(i) = int(tmp(tmp(i,2),1),kind=4)
    end do
    ind = int(tmp(:,2),kind=4)

  end subroutine sortInd



  subroutine sort_global_cols(A)
    use m_mrgrnk
    type(crsMat_t), intent(inout) :: A
    integer(kind=8) :: i, off, n, maxPerRow
    real(kind=8), allocatable :: val(:)
    integer(kind=8), allocatable :: idx(:,:)
    integer(kind=8) :: j

    maxPerRow = maxval(A%row_offset(2:A%nRows+1)-A%row_offset(1:A%nRows))
    allocate(idx(maxPerRow,2))
    allocate(val(maxPerRow))

    ! we still need to sort the entries in each row by column
    do i = 1, A%nRows
      off = A%row_offset(i)
      n = A%row_offset(i+1) - off
      if( n .gt. 1 ) then
        ! reuse idx for sorting
        idx(1:n,1) = A%global_col_idx( off:off+n-1 )
        val(1:n)   = A%val( off:off+n-1 )

        ! determine order of columns
        call mrgrnk(idx(1:n,1), idx(1:n,2))

        ! copy back sorted arrays
        do j = 1, n
          A%global_col_idx( off+j-1 ) = idx(idx(j,2),1)
          A%val( off+j-1 ) = val(idx(j,2))
        end do
      end if
    end do

  end subroutine sort_global_cols


#ifdef PHIST_HAVE_PARMETIS
  !> rearranges the entries in crsMat_t%val and crsMat_t%global_col_idx in a given order
  !! \param crsMat matrix to reorder
  !! \param new_ind array of size crsMat%nEntries with new position of each block
  pure subroutine reord_val_global_col(crsMat, new_ind)
    type(crsMat_t),intent(in out) :: crsMat
    integer(kind=8),intent(in out) :: new_ind(1:crsMat%nEntries)

    integer(kind=8) :: j,l
    integer(kind=8) :: buff_col_ind
    integer(kind=8) :: buff_new_ind
    real(kind=8) :: buff_value

    ! now sort the data in place (we have already the future position,
    ! so this can be done by O(n) swaps)
    do j=1,crsMat%nEntries, 1
      do while( new_ind(j) .ne. j )
        l = new_ind(j)
        ! swap j and l

        buff_value = crsMat%val(l)
        buff_col_ind = crsMat%global_col_idx(l)
        buff_new_ind = new_ind(l)

        crsMat%val(l) = crsMat%val(j)
        crsMat%global_col_idx(l) = crsMat%global_col_idx(j)
        new_ind(l) = new_ind(j)

        crsMat%val(j) = buff_value
        crsMat%global_col_idx(j) = buff_col_ind
        new_ind(j) = buff_new_ind

      end do
    end do

  end subroutine reord_val_global_col
#endif



  !> rearrange the elements of a matrix stored in crs format in a
  !! way that the elements are grouped by remote process (e.g. when performing
  !! a matrix vector multiplication, the elements that are multiplied with a
  !! subvector from a remote process are stored consecutivly in memory)
  subroutine sort_rows_local_nonlocal(mat)
    use m_mrgrnk
    type(CrsMat_t), intent(inout) :: mat

    integer :: i, n
    integer(kind=8) :: j, off, maxPerRow
    integer(kind=8) :: firstrow, lastrow
    integer(kind=8), allocatable :: idx(:,:)
    real(kind=8), allocatable :: val(:)


    firstrow = mat%row_map%distrib(mat%row_map%me)
    lastrow  = mat%row_map%distrib(mat%row_map%me+1)-1

    ! distinguish between local elements and buffer indices
    allocate(mat%nonlocal_offset(mat%nRows))
    mat%nonlocal_offset(:) = mat%row_offset(1:mat%nRows)
!!$omp parallel do schedule(static)
    do i = 1, mat%nRows, 1
      do j = mat%row_offset(i), mat%row_offset(i+1)-1, 1
        if( mat%global_col_idx(j) .ge. firstrow .and. &
          & mat%global_col_idx(j) .le. lastrow        ) then
          ! local element, substract offset of this process
          mat%global_col_idx(j) = (mat%global_col_idx(j)-firstrow+1)
          ! increase nonlocal offset
          mat%nonlocal_offset(i) = mat%nonlocal_offset(i) + 1
        else
          ! nonlocal element, search its position in the comm buffer
          ! add nRows to sort them behind local elements
          mat%global_col_idx(j) = search_recvBuffIndex(mat%global_col_idx(j)) + mat%nRows
        end if
      end do
    end do

    maxPerRow = maxval(mat%row_offset(2:mat%nRows+1)-mat%row_offset(1:mat%nRows))
    allocate(idx(maxPerRow,2))
    allocate(val(maxPerRow))
    allocate(mat%col_idx(mat%nEntries))
    ! sort all rows, try to respect numa
!!$omp parallel do schedule(static) private(idx,val)
    do i = 1, mat%nRows
      off = mat%row_offset(i)
      n = int(mat%row_offset(i+1) - off,kind=4)
      if( n .gt. 0 ) then
        ! reuse idx for sorting
        idx(1:n,1) = mat%global_col_idx( off:off+n-1 )
        val(1:n)   = mat%val( off:off+n-1 )

        ! determine order of columns
        call mrgrnk(idx(1:n,1), idx(1:n,2))

        ! copy back sorted arrays
        do j = 1, n
          mat%col_idx( off+j-1 ) = int(idx(idx(j,2),1),kind=4)
          mat%val( off+j-1 ) = val(idx(j,2))
        end do

      end if
      ! subtract mat%nRows from buffer indices, was only added for sorting!
      do j = mat%nonlocal_offset(i), mat%row_offset(i+1)-1, 1
        mat%col_idx(j) = mat%col_idx(j) - mat%nRows
      end do
    end do
    deallocate(mat%global_col_idx)

  contains

  function search_recvBuffIndex(col) result(idx)
    integer(kind=8), intent(in) :: col
    integer(kind=8) :: idx
    integer(kind=8) :: a, b, n, n_max


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

  end subroutine sort_rows_local_nonlocal


#ifdef PHIST_HAVE_PARMETIS

  !> repartitions a matrix stored in crs form (crsMat_t) in order
  !! to reduce the mpi-communication needed for matrix-vector multiplication
  !! \param crsMat local part of the distributed matrix (gets modified)
  !! \param outlev output verbosity level
  !! \param ierr returns value != 0 on error
  subroutine repartcrs(crsMat,outlev,ierr)
    use mpi
    type(crsMat_t) :: crsMat
    integer,intent(in) :: outlev
    integer,intent(out) :: ierr

    integer(kind=8),allocatable :: rowSendProc(:)
    integer,allocatable :: rowRecvProc(:)
    integer(kind=8),allocatable :: rowSendProcGlob(:)
    integer :: nd_new, jProc, jProc_
    integer :: nSendProcs, nRecvProcs
    integer(kind=8) :: j, k, l, i, row
    integer(kind=8),allocatable :: sendBuffInd(:)
    integer(kind=8),allocatable :: sendIds(:),procSendCount(:)
    integer(kind=8),allocatable :: sendIdsInd(:)
    integer,allocatable :: recvRequest(:,:)
    integer,allocatable :: recvStatus(:,:,:)
    integer(kind=8),allocatable :: sendBuff_row_blk(:)
    integer(kind=8),allocatable :: buff_ind(:)
    integer(kind=8),allocatable :: col_ind_new(:)
    integer(kind=8),allocatable :: row_blk_new(:)
    real(kind=8),allocatable :: value_new(:)
    integer(kind=8),allocatable :: glob_row_permut(:)
    integer,allocatable :: offsets(:)
    integer :: ierr_, localErr


    ierr = 0


    allocate(rowSendProc(crsMat%nRows))
    ! don't call parmetis when we only have very vew elements per process, it'll segfault otherwise somehow...
    if( (crsMat%row_map%distrib(crsMat%row_map%nProcs)-1) / crsMat%row_map%nProcs .lt. 10 ) then
      rowSendProc = crsMat%row_map%me
      if( outlev .ge. 1 .and. crsMat%row_map%me .eq. 0 ) then
        write(*,*) 'WARNING: disabling parmetis as we have too few elements per node!'
      end if
    else
      ! parmetis call encapsulated, we only need the
      ! data where each row goes
      call calculateNewPartition(crsMat,rowSendProc)
    end if
!write(*,*) 'rowSendProc', rowSendProc
    localErr = 0
    if( any(rowSendProc .lt. 0 .or. rowSendProc .ge. crsMat%row_map%nProcs) ) then
      localErr = 1
    end if
    call mpi_allreduce(localErr, ierr, 1, MPI_INTEGER, MPI_SUM, crsMat%row_map%comm, ierr_)
    if( ierr .ne. 0 ) then
      if( outlev .ge. 1 .and. crsMat%row_map%me .eq. 0 ) then
        write(*,*) 'ERROR: parmetis failure!'
      end if
      return
    end if

    ! determine recvRowProc vector
    allocate(rowSendProcGlob(crsMat%row_map%distrib(crsMat%row_map%nProcs)-1),&
             offsets(crsMat%row_map%nProcs), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
#warning "cannot allgatherv more than max-int32 elements!"
    offsets = crsMat%row_map%distrib(0:crsMat%row_map%nProcs-1)-1
!write(*,*) 'offsets', offsets
    call mpi_allgatherv(rowSendProc,crsMat%nRows,MPI_INTEGER8,&
      & rowSendProcGlob,crsMat%row_map%nlocal,offsets,MPI_INTEGER8, &
      & crsMat%row_map%comm, ierr)
      if (ierr/=0) then
        return
      end if
!if( crsMat%row_map%me .eq. 0 ) then
  !write(*,*) 'rowSendProcGlob', rowSendProcGlob
!end if
    nd_new = 0
    do i = 1, crsMat%row_map%distrib(crsMat%row_map%nProcs)-1, 1
      if( rowSendProcGlob(i) .eq. crsMat%row_map%me ) then
        nd_new = nd_new + 1
      end if
    end do

    allocate(rowRecvProc(nd_new), &
             crsMat%global_row_idx(nd_new), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    nRecvProcs = 1
    if( nd_new .eq. 0 ) nRecvProcs = 0
    nd_new = 0
    do jProc = 0, crsMat%row_map%nProcs-1, 1
      do i = crsMat%row_map%distrib(jProc), crsMat%row_map%distrib(jProc+1)-1, 1
        if( rowSendProcGlob(i) .eq. crsMat%row_map%me ) then
          nd_new = nd_new + 1
          rowRecvProc(nd_new) = jProc
          if( nd_new .gt. 1 ) then
            if( rowRecvProc(nd_new-1) .ne. jProc ) then
              nRecvProcs = nRecvProcs + 1
            end if
          end if
          crsMat%global_row_idx(nd_new) = i
        end if
      end do
    end do
    ! copy to map
    allocate(crsMat%row_map%global_idx(size(crsMat%global_row_idx)))
    crsMat%row_map%global_idx = crsMat%global_row_idx

    crsMat%row_map%distrib=0
    do i = 1, size(rowSendProcGlob), 1
      crsMat%row_map%distrib(rowSendProcGlob(i)+1) = &
        crsMat%row_map%distrib(rowSendProcGlob(i)+1) + 1
    end do
    crsMat%row_map%distrib(0) = 1
    do i = 1, crsMat%row_map%nProcs, 1
      crsMat%row_map%distrib(i) = crsMat%row_map%distrib(i-1) + crsMat%row_map%distrib(i)
    end do
    crsMat%row_map%nlocal = int(crsMat%row_map%distrib(1:crsMat%row_map%nProcs) - crsMat%row_map%distrib(0:crsMat%row_map%nProcs-1),kind=4)
    allocate(glob_row_permut(crsMat%row_map%distrib(crsMat%row_map%nProcs)-1), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    do i = 1, size(glob_row_permut), 1
      glob_row_permut( i ) = crsMat%row_map%distrib( rowSendProcGlob(i) )
      crsMat%row_map%distrib( rowSendProcGlob(i) ) = &
        & crsMat%row_map%distrib( rowSendProcGlob(i) ) + 1
    end do
    ! restore row_blk_offset
    do i = crsMat%row_map%nProcs,1,-1
      crsMat%row_map%distrib(i) = crsMat%row_map%distrib(i-1)
    end do
    crsMat%row_map%distrib(0) = 1
    deallocate(rowSendProcGlob)

!write(*,*) 'global_col_idx', crsMat%global_col_idx
!write(*,*) 'glob_row_permut', glob_row_permut
    ! change the col_ind for the new partition
    do j = 1, crsMat%nEntries, 1
      crsMat%global_col_idx(j) = glob_row_permut( crsMat%global_col_idx(j) )
    end do
    deallocate(glob_row_permut)
    call sort_global_cols(crsMat)


    ! we now have the rowRecvProc and the rowSendProc vec.
    ! first start the receive for the number of elemens per row
    allocate(recvRequest(nRecvProcs,2))
    recvRequest=MPI_REQUEST_NULL
    allocate(row_blk_new(nd_new+1))
    if( nd_new .gt. 0 ) jProc = rowRecvProc(1)
    k = 1
    l = 1
    row_blk_new = 0
    do i = 2, nd_new+1, 1
      jProc_ = -1
      if( i .le. nd_new ) jProc_ = rowRecvProc(i)
      if( jProc .ne. jProc_ ) then
        call mpi_irecv(row_blk_new(k+1),int(i-k),MPI_INTEGER8, &
          &           jProc,0,crsMat%row_map%comm, recvRequest(l,1), ierr)
        k = i
        l = l + 1
        jProc = jProc_
      end if
    end do


    ! now create the send buffer and send the number of elements per row
    allocate(procSendCount(0:crsMat%row_map%nProcs))
    procSendCount = 0
    nSendProcs = 0
    do i = 1, crsMat%nRows, 1
      if( procSendCount(rowSendProc(i)+1) .eq. 0 ) then
        nSendProcs = nSendProcs + 1
      end if
      procSendCount(rowSendProc(i)+1) = procSendCount(rowSendProc(i)+1) + 1
    end do
    allocate(sendBuffInd(nSendProcs+1))
    allocate(sendIds(nSendProcs))
    allocate(sendIdsInd(0:crsMat%row_map%nProcs-1))
    jProc = 0
    sendBuffInd(1) = 0
    sendIdsInd=-1
    do i = 1, nSendProcs, 1
      do while( procSendCount(jProc+1) .eq. 0 )
        jProc = jProc + 1
      end do
      sendBuffInd(i+1) = procSendCount(jProc+1)
      sendIds(i) = jProc
      sendIdsInd(jProc) = i
      jProc=jProc+1
    end do
    sendBuffInd(1) = 1
    do i = 1, nSendProcs, 1
      sendBuffInd(i+1) = sendBuffInd(i)+sendBuffInd(i+1)
    end do
    do i = 1, crsMat%row_map%nProcs, 1
      procSendCount(i) = procSendCount(i)+procSendCount(i-1)
    end do
    allocate(sendBuff_row_blk(sendBuffInd(nSendProcs+1)), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    sendBuff_row_blk=0
    do i = 1, crsMat%nRows, 1
      rowSendProc(i) = sendIdsInd(rowSendProc(i))
      procSendCount(sendIds(rowSendProc(i))) = &
        & procSendCount(sendIds(rowSendProc(i)))+1
      sendBuff_row_blk( procSendCount(sendIds(rowSendProc(i))) + 1) = &
        & crsMat%row_offset(i+1) - crsMat%row_offset(i)
    end do
    deallocate(sendIdsInd)
    do i = crsMat%row_map%nProcs, 1, -1
      procSendCount(i) = procSendCount(i-1)
    end do
    procSendCount(0) = 0
!write(*,*) 'row_offset', crsMat%row_offset
!write(*,*) 'sendBuff_row_blk', sendBuff_row_blk
    do i = 1, nSendProcs, 1
      call mpi_send(sendBuff_row_blk(sendBuffInd(i)+1), &
        &     int(sendBuffInd(i+1)-sendBuffInd(i)), MPI_INTEGER8, &
        &     int(sendIds(i)), 0, crsMat%row_map%comm, ierr)
    end do
    sendBuff_row_blk(1) = 1
    do i = 1, sendBuffInd(nSendProcs+1)-1, 1
      sendBuff_row_blk(i+1) = sendBuff_row_blk(i) + sendBuff_row_blk(i+1)
    end do

    ! finish receive and allocate buffers for col_ind receive
    allocate(recvStatus(MPI_STATUS_SIZE,nRecvProcs,2))
    call mpi_waitall(nRecvProcs,recvRequest(:,1),recvStatus(:,:,1),ierr)
    row_blk_new(1) = 1
    do i = 1, nd_new, 1
      row_blk_new(i+1) = row_blk_new(i) + row_blk_new(i+1)
    end do
!write(*,*) 'row_blk_new', row_blk_new
    allocate(col_ind_new(row_blk_new(nd_new+1)-1), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    if( nd_new .gt. 0 ) jProc = rowRecvProc(1)
    k = 1
    l = 1
    do i = 2, nd_new+1, 1
      jProc_ = -1
      if( i .le. nd_new ) jProc_ = rowRecvProc(i)
      if( jProc .ne. jProc_ ) then
        call mpi_irecv(col_ind_new(row_blk_new(k)),&
          &           int(row_blk_new(i)-row_blk_new(k)),MPI_INTEGER8, &
          &           jProc,1,crsMat%row_map%comm, recvRequest(l,1), ierr)
        k = i
        l = l + 1
        jProc = jProc_
      end if
    end do




    ! reorder local matrix entries, so they can be sent block-wise
    allocate(buff_ind(crsMat%nEntries), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    do i = 1, crsMat%nRows, 1
      procSendCount(sendIds(rowSendProc(i))) = &
        & procSendCount(sendIds(rowSendProc(i))) + 1
      row = procSendCount(sendIds(rowSendProc(i)))
      do j = crsMat%row_offset(i), crsMat%row_offset(i+1)-1, 1
        buff_ind(j) = sendBuff_row_blk( row )
        sendBuff_row_blk(row) = sendBuff_row_blk(row) + 1
      end do
    end do
    do i = crsMat%nRows+1, 2, -1
      sendBuff_row_blk(i) = sendBuff_row_blk(i-1)      
    end do
    sendBuff_row_blk(1) = 1
    call reord_val_global_col(crsMat,buff_ind)
    deallocate(buff_ind,procSendCount,crsMat%row_offset)

!write(*,*) 'global_col_idx', crsMat%global_col_idx
    ! send col ind
    do i = 1, nSendProcs, 1
      call mpi_send(crsMat%global_col_idx(sendBuff_row_blk(sendBuffInd(i))), &
        &  int(sendBuff_row_blk(sendBuffInd(i+1))-sendBuff_row_blk(sendBuffInd(i))), &
        &  MPI_INTEGER8,int(sendIds(i)),1,crsMat%row_map%comm, ierr)
    end do
    deallocate(crsMat%global_col_idx) ! not needed any more


    ! finish col_ind receive and create buffer for value
    call mpi_waitall(nRecvProcs,recvRequest(:,1),recvStatus(:,:,1),ierr)
    allocate(value_new(row_blk_new(nd_new+1)-1), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    if( nd_new .gt. 0 ) jProc = rowRecvProc(1)
    k = 1
    l = 1
    do i = 2, nd_new+1, 1
      jProc_ = -1
      if( i .le. nd_new ) jProc_ = rowRecvProc(i)
      if( jProc .ne. jProc_ ) then
        call mpi_irecv(value_new(row_blk_new(k)),&
          &           int(row_blk_new(i)-row_blk_new(k)),&
          &           MPI_DOUBLE_PRECISION, &
          &           jProc,2,crsMat%row_map%comm, recvRequest(l,1), ierr)
        k = i
        l = l + 1
        jProc = jProc_
      end if
    end do


    ! send value
    do i = 1, nSendProcs, 1
      call mpi_send(crsMat%val(sendBuff_row_blk(sendBuffInd(i))), &
        &  int(sendBuff_row_blk(sendBuffInd(i+1))-sendBuff_row_blk(sendBuffInd(i))), &
        &  MPI_DOUBLE_PRECISION,int(sendIds(i)),2,crsMat%row_map%comm, ierr)
    end do
    deallocate(crsMat%val,sendBuff_row_blk) ! not needed any more


    ! finish value receive
    call mpi_waitall(nRecvProcs,recvRequest(:,1),recvStatus(:,:,1),ierr)
    deallocate(sendBuffInd, sendIds, rowRecvProc, rowSendProc,recvRequest,recvStatus)




    ! copy resulting data back

    ! try to respect NUMA
    crsMat%nRows = nd_new
    crsMat%nEntries = row_blk_new(nd_new+1)-1
    allocate(crsMat%row_offset(crsMat%nRows+1), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
!$omp parallel do schedule(static)
    do i = 1, crsMat%nRows+1, 1
      crsMat%row_offset(i) = row_blk_new(i)
    end do
    deallocate(row_blk_new)

    allocate(crsMat%global_col_idx(crsMat%nEntries), &
             crsMat%val(crsMat%nEntries), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if             
!$omp parallel do schedule(static)
    do i = 1, crsMat%nRows, 1
      do j = crsMat%row_offset(i), crsMat%row_offset(i+1)-1, 1
        crsMat%global_col_idx(j) = col_ind_new(j)
        crsMat%val(j) = value_new(j)
      end do
    end do
    deallocate(col_ind_new)
    deallocate(value_new)


  contains


    subroutine calculateNewPartition(crsMat, rowBlkTargetProc)
      use iso_c_binding
      type(crsMat_t),intent(in) :: crsMat
      integer(kind=C_INT64_T) :: rowBlkTargetProc(crsMat%nRows)


      ! parameter for ParMETIS_V3_PartKway
      integer(kind=C_INT64_T),allocatable :: vtxdist(:) ! = crsMat%row_map%distrib+1
      integer(kind=C_INT64_T),allocatable :: xadj(:) ! = crsMat%row_offset without counting diagonal entries
      integer(kind=C_INT64_T),allocatable :: adjncy(:) ! = crsMat%global_col_idx without diagonal (connection to itself)
      integer(kind=C_INT64_T),allocatable :: vwgt(:) ! has size(crsMat%nRows), number of non-zero entries per row block
      integer(kind=C_INT64_T),allocatable :: adjwgt(:) ! edge weights, for unsymmetric matrices 1 or 2
      integer(kind=C_INT64_T),parameter :: wgtflag = 2 ! vertex and edge weights
      integer(kind=C_INT64_T),parameter :: numflag = 1 ! fortran numbering style
      integer(kind=C_INT64_T),parameter :: ncon = 1 ! only one constraint, one weight per vertex
      integer(kind=C_INT64_T) :: nparts ! = crsMat%row_map%nProcs
      real(kind=C_FLOAT),allocatable :: tpwgts(:) ! =1/nparts, has size (1:nparts)
      real(kind=C_FLOAT),parameter :: ubvec = 1.01 ! vertex weight imbalance tolerance
      real(kind=C_FLOAT),parameter :: ubvec_refine = 1.01 ! vertex weight imbalance tolerance
      integer(kind=C_INT64_T) :: options(0:3) = 0 ! default options
      integer(kind=C_INT64_T) :: edgecut_refined, edgecut ! output parameter, global number of cutted edges
      integer(kind=C_INT64_T) :: edgecut_before


      integer :: i, ierr
      integer(kind=8) :: j, k
      integer(kind=8) :: offset

      offset = crsMat%row_map%distrib(crsMat%row_map%me)-1
      allocate(vtxdist(0:crsMat%row_map%nProcs))
      vtxdist=crsMat%row_map%distrib
      ! count non diagonal elements
      allocate(xadj(crsMat%nRows+1), &
               vwgt(crsMat%nRows), stat=ierr)
      if (ierr/=0) then
        ierr=-44
        return
      end if
      j = 1
      xadj(1) = 1
      do i = 1, crsMat%nRows, 1
        do k = crsMat%row_offset(i), crsMat%row_offset(i+1)-1, 1
          if( crsMat%global_col_idx(k) .ne. i + offset ) then
            j = j + 1
          end if
        end do
        xadj(i+1) = j
        vwgt(i) = crsMat%row_offset(i+1)-crsMat%row_offset(i)
      end do
      allocate(adjncy(xadj(crsMat%nRows+1)-1), stat=ierr)
      if (ierr/=0) then
        ierr=-44
        return
      end if
      !nullify(adjwgt)
      j = 1
      edgecut = 0
      do i = 1, crsMat%nRows, 1
        do k = crsMat%row_offset(i), crsMat%row_offset(i+1)-1, 1
          if( crsMat%global_col_idx(k) .ne. i + offset ) then
            adjncy(j) = crsMat%global_col_idx(k)
            j = j + 1
            if( crsMat%global_col_idx(k) .le. offset .or. &
              & crsMat%global_col_idx(k) .gt. offset+crsMat%nRows ) then
              edgecut=edgecut+1
            end if
          end if
        end do
      end do

      call mpi_allreduce(edgecut,edgecut_before,1,MPI_INTEGER8,MPI_SUM,crsMat%row_map%comm,ierr)

      nparts = crsMat%row_map%nProcs
      allocate(tpwgts(1:nparts))
      tpwgts = 1.0/nparts
      rowBlkTargetProc=0

      ! make symmetric graph from non-symmetric matrix pattern
      !call parudgraph(mpi,vtxdist,xadj,adjncy,adjwgt)
      allocate(adjwgt(xadj(crsMat%nRows+1)-1),stat=ierr)
      if (ierr/=0) then
        ierr=-44
        return
      end if
      adjwgt = 2
      ! the edge weights are 1 for single edges and 2 for edges in both
      ! directions

!write(*,*) 'vtxdist', vtxdist
!write(*,*) 'xadj', xadj
!write(*,*) 'adjncy', adjncy
!write(*,*) 'vwgt', vwgt
!write(*,*) 'adjwgt', adjwgt
!write(*,*) 'wgtflag', wgtflag
!write(*,*) 'numflag', numflag
!write(*,*) 'ncon', ncon
!write(*,*) 'nparts', nparts
!write(*,*) 'tpwgts', tpwgts
!write(*,*) 'ubvec', ubvec
!write(*,*) 'options', options
!write(*,*) 'edgecut', edgecut
!write(*,*) 'rowBlkTargetProc', rowBlkTargetProc
      call ParMETIS_V3_PartKway_f(vtxdist, xadj, adjncy, vwgt, adjwgt, &
        &                         wgtflag, numflag, ncon, nparts, tpwgts, ubvec, &
        &                         options, edgecut, rowBlkTargetProc, crsMat%row_map%comm)
!write(*,*) 'rowBlkTargetProc', rowBlkTargetProc

      ! doesn't work with ParMetis-3.1.1, but other versions fail for empty
      ! processes
      if( outlev .ge. 3 .and. crsMat%row_map%me .eq. 0 ) then
        write(*,*) 'repartitioning: total number of edges cut:'
        write(*,*) '                before:',edgecut_before
        write(*,*) '      after first step:',edgecut
      end if

      if( edgecut .gt. 0 ) then
        options(0) = 1
        options(1) = 0
        options(2) = 15
        options(3) = 2

        call ParMETIS_V3_RefineKway_f(vtxdist, xadj, adjncy, vwgt, adjwgt, &
          &                           wgtflag, numflag, ncon, nparts, tpwgts, ubvec_refine, &
          &                           options, edgecut_refined, rowBlkTargetProc, crsMat%row_map%comm)

        if( outlev .ge. 3 .and. crsMat%row_map%me .eq. 0 ) then
          write(*,*) '     after second step:',edgecut_refined
        end if
      end if

      deallocate(xadj,adjncy,vwgt,tpwgts,vtxdist,adjwgt)

      rowBlkTargetProc = rowBlkTargetProc - 1

    end subroutine calculateNewPartition
  end subroutine repartcrs



#endif /* PHIST_HAVE_PARMETIS */

#ifdef PHIST_HAVE_COLPACK

  subroutine colorcrs(crsMat,dist,outlev,ierr)
      use, intrinsic :: iso_c_binding
    use mpi
    implicit none
    type(crsMat_t),target :: crsMat
    ! distance (1 or 2)
    integer, intent(in) :: dist
    integer,intent(in) :: outlev
    integer,intent(out) :: ierr
    
    ! local variables
    integer :: i
    integer(kind=8) :: j
    integer(c_int), target, allocatable :: colors(:)
    integer :: ncolors
    integer, allocatable, dimension(:) :: colorCount
    
    TYPE(c_ptr) :: c_rp, c_nlp, c_cidx, c_colors
    
    interface
    
      subroutine colpack(nloc,row_ptr,nonlocal_ptr,column_idx, &
              ncolors, colors, dist, idxBase,ierr) &
        bind(C, name='colpack_v1_0_8')

      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(c_int), value :: nloc
      type(c_ptr), value :: row_ptr, nonlocal_ptr, column_idx
      integer(c_int), intent(out) :: ncolors      
      type(c_ptr), value :: colors
      integer(c_int), value :: dist
      integer(c_int), value :: idxBase
      integer(c_int), intent(inout) :: ierr

      end subroutine colpack

    end interface

    allocate(colors(crsMat%nRows))

    c_rp=c_loc(crsMat%row_offset(1))
    c_nlp=c_loc(crsMat%nonlocal_offset(1))
    c_cidx=c_loc(crsMat%col_idx(1))
    c_colors=c_loc(colors(1))
        
      call colpack(crsMat%nRows, c_rp,c_nlp,c_cidx,ncolors,c_colors, dist,1,ierr)
      
      if (ierr/=0) then
        write(*,*) 'error code ',ierr,' returned from colpack'
        return
      end if
 
    ierr=0
    
    ! adjust colors to 1-based indexing
    colors(1:crsMat%nRows)=colors(1:crsMat%nRows)+1

    ! count number of members in each color set
    allocate(colorCount(ncolors))
    colorCount(1:ncolors)=0
    do i=1,crsMat%nRows
      colorCount(colors(i))=colorCount(colors(i))+1
    end do
    ! allocate the color arrays in the map of A and fill them.
    ! They are deleted when the map is deleted.
    crsMat%row_map%coloringType=2
    crsMat%row_map%nColors=ncolors
    allocate(crsMat%row_map%color_offset(ncolors+1))
    allocate(crsMat%row_map%color_idx(crsMat%nRows))
    
    crsMat%row_map%color_offset(1)=1
    do i=1,ncolors
      crsMat%row_map%color_offset(i+1) = &
      crsMat%row_map%color_offset(i)   + colorCount(i)
    end do

    colorCount(1:ncolors)=0
    do i=1,crsMat%nRows
      j=crsMat%row_map%color_offset(colors(i))+colorCount(colors(i))
      colorCount(colors(i))=colorCount(colors(i))+1
      crsMat%row_map%color_idx(j)=i
    end do

    deallocate(colorCount)
#if defined(TESTING) && PHIST_OUTLEV >= 4
#define DEBUG_COLPACK 1
#endif
#if DEBUG_COLPACK
! for n>=10,000 the file gets big and the format used below won't fit
if (crsMat%nRows<10000) then
        open(unit=42,file='test_coloring.m',status='replace')
        write(42,*) 'n=',crsMat%nRows
        write(42,*) 'Adat=[...'
        do i=1,crsMat%nRows
          do j=crsMat%row_offset(i),crsMat%nonlocal_offset(i)-1
            write(42,"(I4,I4,E24.8)") i,crsMat%col_idx(j),crsMat%val(j)
          end do
        end do
        write(42,*) '];'
        write(42,*) 'A=sparse(Adat(:,1),Adat(:,2),Adat(:,3),n,n);'
        write(42,*) 'ncolors=',ncolors,';'
        write(42,*) 'colors=[...'
        do i=1,crsMat%nRows
          write(42,*) colors(i)
        end do
        write(42,*) '];'
        write(42,*) 'color_ptr=[...'
                do i=1,crsMat%row_map%nColors+1
          write(42,*) crsMat%row_map%color_offset(i)
        end do
        write(42,*) '];'
        write(42,*) 'color_idx=[...'
                do i=1,crsMat%nRows
          write(42,*) crsMat%row_map%color_idx(i)
        end do
        write(42,*) '];'

!       permute local CRS matrix rows and cols, this will mess up MPI parallel runs and
!       result in permuted vectors, but just for single-node performance we do it right now
!call permute_local_matrix(crsMat)

        write(42,*) 'Adat2=[...'
        do i=1,crsMat%nRows
          do j=crsMat%row_offset(i),crsMat%nonlocal_offset(i)-1
            write(42,"(I4,I4,E24.8)") i,crsMat%col_idx(j),crsMat%val(j)
          end do
        end do
        write(42,*) '];'
        write(42,*) 'A2=sparse(Adat2(:,1),Adat2(:,2),Adat2(:,3),n,n);'

        flush(42)
        close(42)
end if
#else
!       permute local CRS matrix rows and cols, this will mess up MPI parallel runs and
!       result in permuted vectors, but just for single-node performance we do it right now
!call permute_local_matrix(crsMat)
#endif

    deallocate(colors)

    end subroutine colorcrs

subroutine permute_local_matrix(A)

implicit none
TYPE(crsMat_t) :: A
integer, allocatable, dimension(:) :: row_ptr
integer, allocatable, dimension(:) :: nonlocal_ptr
integer, allocatable, dimension(:) :: col_idx
real(kind=8), allocatable, dimension(:) :: val
integer, allocatable, dimension(:) :: inv_perm

integer :: ic,jc,i,k,l,ii,len_i, len_nl_i
integer(kind=8) :: j,jj

allocate(inv_perm(A%nRows))
allocate(row_ptr(A%nRows+1))
allocate(nonlocal_ptr(A%nRows))
allocate(val(A%nEntries))
allocate(col_idx(A%nEntries))

do i=1,A%nRows
  inv_perm(A%row_map%color_idx(i))=i
end do

#ifdef DEBUG_COLPACK
write(42,*) 'iperm=[...'
do i=1,A%nRows
  write(42,*) inv_perm(i)
end do
write(42,*) '];'
#endif

#ifdef TESTING
! quick test if the pointer is correct
if (A%row_map%color_offset(A%row_map%nColors+1)-1/=A%nRows) then
  stop 'color_ptr incorrect'
end if
do i=1,A%nRows
  if (A%row_map%color_idx(i)<=0 .or. &
      A%row_map%color_idx(i)>A%nRows) then
    stop 'color_idx incorrect'
  end if
end do
#endif

row_ptr(1)=1
do ic=1,A%row_map%nColors
  ! new row index ii: colors contiguous
  do ii=A%row_map%color_offset(ic),A%row_map%color_offset(ic+1)-1
    ! old row index i: where to find row ii in old matrix
    i=A%row_map%color_idx(ii)
    len_i= A%nonlocal_offset(i)-A%row_offset(i)
    len_nl_i=A%row_offset(i+1)-A%nonlocal_offset(i)
    row_ptr(ii+1)=row_ptr(ii) + len_i
    nonlocal_ptr(ii)=row_ptr(ii+1)+len_nl_i
    do j=A%row_offset(i),A%nonlocal_offset(i)-1
      jj=row_ptr(ii)+ (j-A%row_offset(i))
      col_idx(jj)=inv_perm(A%col_idx(j))
      val(jj)=A%val(j)
    end do
    do j=A%nonlocal_offset(i),A%row_offset(i+1)-1
      jj=row_ptr(ii)+ (j-A%row_offset(i))
      col_idx(jj)=A%col_idx(j)
      val(jj)=A%val(j)
    end do
  end do
end do

! permute comm buffer indices
do i=1,A%comm_buff%nSendProcs, 1
  k = A%comm_buff%sendInd(i)
  l = A%comm_buff%sendInd(i+1)
  do j = k, l-1, 1
    !send elt j comes from x(sendRowBlkInd(j)
    A%comm_buff%sendRowBlkInd(j) = &
    A%row_map%color_idx(A%comm_buff%sendRowBlkInd(j))
  end do
end do

do i=1,A%comm_buff%nRecvProcs, 1
  k = A%comm_buff%recvInd(i)
  l = A%comm_buff%recvInd(i+1)
  do j = k, l-1, 1
    !send elt j comes from x(sendRowBlkInd(j)
    A%comm_buff%recvRowBlkInd(j) = &
    A%row_map%color_idx(A%comm_buff%recvRowBlkInd(j))
  end do
end do
!write(*,*) 'dat0=[...'
!do i=1,A%nRows
!  do j=A%row_offset(i),A%nonlocal_offset(i)
!    write(*,*) i,A%col_idx(j),A%val(j)
!  end do
!end do
!write(*,*) '];'
A%row_offset(1:A%nRows+1) = row_ptr(1:A%nRows+1)
A%nonlocal_offset(1:A%nRows) = nonlocal_ptr(1:A%nrows)
A%col_idx(1:A%nEntries) = col_idx(1:A%nEntries)
A%val(1:A%nEntries) = val(1:A%nEntries)

deallocate(inv_perm)
deallocate(row_ptr)
deallocate(nonlocal_ptr)
deallocate(col_idx)
deallocate(val)

! TODO: permute comm buffer indices
!write(*,*) 'dat1=[...'
!do i=1,A%nRows
!  do j=A%row_offset(i),A%nonlocal_offset(i)
!    write(*,*) i,A%col_idx(j),A%val(j)
!  end do
!end do
!write(*,*) '];'
end subroutine permute_local_matrix

#endif

  !==================================================================================
  !> multiply crsmat with mvec
  subroutine crsmat_times_mvec(alpha, A, shifts, x, beta, y, ierr)
    use, intrinsic :: iso_c_binding, only: C_INT
    use mpi
    !--------------------------------------------------------------------------------
    real(kind=8),   intent(in)    :: alpha
    type(CrsMat_t), intent(inout) :: A
    real(kind=8),   intent(in)    :: shifts(*)
    type(MVec_t),   intent(in)    :: x
    real(kind=8),   intent(in)    :: beta
    type(MVec_t),   intent(inout) :: y
    integer(C_INT), intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    integer :: nvec, ldx, ldy, recvBuffSize, sendBuffSize
    logical :: strided_x, strided_y, strided
    logical :: y_is_aligned16, handled
    integer :: i, k, l
    !--------------------------------------------------------------------------------

    ! check arguments
    ierr = 0
    if( y%jmax-y%jmin .ne. x%jmax-x%jmin ) then
      ierr = -1
      return
    end if
#ifdef TESTING
    if( .not. map_compatible_map(A%row_map, x%map) .or. &
      & .not. map_compatible_map(A%row_map, y%map)      ) then
      ierr = -1
      return
    end if
#endif

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

#if defined(TESTING) && PHIST_OUTLEV >= 4
    write(*,*) 'spMVM with nvec =',nvec,', ldx =',ldx,', ldy =',ldy,', y_mem_aligned16 =', y_is_aligned16, 'on proc', A%row_map%me
    flush(6)
#endif


    ! exchange necessary x data => A%comm_buff%recvData
    call Dspmv_import(A,x,ierr)

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

    !--------------------------------------------------------------------------------
  end subroutine crsmat_times_mvec


  !==================================================================================
  !> read MatrixMarket file
  subroutine phist_DcrsMat_read_mm(A_ptr, comm_ptr, filename_len, filename_ptr, ierr) &
    & bind(C,name='phist_DcrsMat_read_mm_f') ! circumvent bug in opari (openmp instrumentalization)
    use, intrinsic :: iso_c_binding
    use env_module, only: newunit
    use mpi
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: A_ptr
    type(C_PTR),        value  :: comm_ptr
    integer, pointer :: comm
    integer(C_INT),     value       :: filename_len
    character(C_CHAR),  intent(in)  :: filename_ptr(filename_len)
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    character(len=filename_len) :: filename
    !--------------------------------------------------------------------------------
    integer :: funit
    character(len=100) :: line
    integer(kind=8), allocatable :: idx(:,:)
    real(kind=8), allocatable :: val(:)
    integer(kind=8) :: i, j, k, globalRows, globalCols, globalEntries, nRows
    !integer(kind=8) :: i_
    logical :: symmetric
    integer(kind=8) :: tmp_idx(2,2)
    real(kind=8) :: tmp_val(2)
    integer :: newEntries
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
    if( trim(line) .eq. '%%MatrixMarket matrix coordinate real general' ) then
      symmetric = .false.
    elseif( trim(line) .eq. '%%MatrixMarket matrix coordinate real symmetric' ) then
      symmetric = .true.
    else
      write(*,*) 'unsupported format'
      flush(6)
      ierr = -99
      return
    end if

    ! read second line (comment)
    read(funit,'(A)') line
    line = adjustl(line)
    do while( line(1:1) .eq. '%' )
      read(funit,'(A)') line
      line = adjustl(line)
    end do

    allocate(A)

    ! now read the dimensions
    read(line,*) globalRows, globalCols, nRows
    if( symmetric ) then
      globalEntries = nRows * 2 - min(globalRows,globalCols)
    else
      globalEntries = nRows
    end if
    write(*,*) 'CrsMat:', globalRows, globalCols, globalEntries
    flush(6)
    call c_f_pointer(comm_ptr, comm)
    call map_setup(A%row_map, comm, globalRows, ierr)
    if( ierr .ne. 0 ) return

    A%nRows = A%row_map%nlocal(A%row_map%me)
    A%nCols = A%nRows

    ! allocate temporary buffers
    allocate(idx(globalEntries,2), &
             val(globalEntries),&
             stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if

    ! read data
    j = 0
    do i = 1, nRows, 1
      read(funit,*) tmp_idx(1,1), tmp_idx(1,2), tmp_val(1)
      if( symmetric .and. tmp_idx(1,1) .ne. tmp_idx(1,2) ) then
        tmp_idx(2,1) = tmp_idx(1,2)
        tmp_idx(2,2) = tmp_idx(1,1)
        tmp_val(2)   = tmp_val(1)
        newEntries = 2
      else
        newEntries = 1
      end if

      do k = 1, newEntries, 1
        idx(j+1,:) = tmp_idx(k,:)
        val(j+1) = tmp_val(k)

        if( idx(j+1,1) .ge. A%row_map%distrib(A%row_map%me) .and. &
          & idx(j+1,1) .lt. A%row_map%distrib(A%row_map%me+1)     ) then
          j = j + 1
          ! already subtract offset of this proc
          idx(j,1) = idx(j,1) -  A%row_map%distrib(A%row_map%me) + 1
        end if
      end do

    end do
    A%nEntries = j

    ! close the file
    close(funit)


    ! allocate crs matrix
    allocate(A%row_offset(A%nRows+1),&
             A%global_col_idx(A%nEntries),&
             A%val(A%nEntries),stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if

    ! try to respect NUMA
!$omp parallel do schedule(static)
    do i = 1, A%nRows+1
      A%row_offset(i) = 0
    end do
!$omp barrier

    ! count number of entries per row
    do i = 1, A%nEntries
      A%row_offset( idx(i,1)+1 ) = A%row_offset( idx(i,1)+1 )  +  1
    end do
    A%row_offset(1) = 1
    do i = 1, A%nRows
      A%row_offset( i+1 ) = A%row_offset( i ) + A%row_offset( i+1 )
    end do

    ! try to respect NUMA
!$omp parallel do schedule(static)
    do i = 1, A%nRows
      do j = A%row_offset(i), A%row_offset(i+1)-1, 1
        A%global_col_idx(j) = 0
        A%val(j) = 0._8
      end do
    end do
!$omp barrier

    ! now put the entries into the corresponding rows
    do i = 1, A%nEntries
      A%global_col_idx( A%row_offset(idx(i,1)) ) = idx(i,2)
      A%val    ( A%row_offset(idx(i,1)) ) = val(i)
      A%row_offset( idx(i,1) ) = A%row_offset( idx(i,1) )  +  1
    end do
    ! set row_offset
    do i = A%nRows, 1, -1
      A%row_offset( i+1 ) = A%row_offset( i )
    end do
    A%row_offset( 1 ) = 1

    call sort_global_cols(A)
#ifdef PHIST_HAVE_PARMETIS
!write(*,*) 'row_offset', A%row_offset
!write(*,*) 'col_idx', A%global_col_idx
!write(*,*) 'val', A%val
    call repartcrs(A,3,ierr)
    if (ierr/=0) then
      return
    end if
!write(*,*) 'row_offset', A%row_offset
!write(*,*) 'col_idx', A%global_col_idx
!write(*,*) 'val', A%val
#endif

! write matrix to mat.mm
!do i_ = 0, A%row_map%nProcs-1

  !call MPI_Barrier(A%row_map%comm, ierr)
  !if( A%row_map%me .ne. i_ ) cycle

  !if( i_ .eq. 0 ) then
    !open(unit   = newunit(funit), file = 'out.mm', action = 'write')
    !write(funit,'(A)') '%%MatrixMarket matrix coordinate real general'
    !write(funit,'(A)') '% generated by SpinChainSZ'
    !write(funit,*) A%row_map%distrib(A%row_map%nProcs)-1, A%row_map%distrib(A%row_map%nProcs)-1, A%nEntries
  !else
    !open(unit   = newunit(funit), file = 'out.mm', action = 'write', position = 'append')
  !end if

  !do i = 1, A%nRows
    !do j = A%row_offset(i), A%row_offset(i+1)-1, 1
      !write(funit,*) A%row_map%distrib(A%row_map%me)+i-1, A%global_col_idx(j), A%val(j)
    !end do
  !end do
  !close(funit)

!end do


    call setup_commBuff(A, A%comm_buff)
    call sort_rows_local_nonlocal(A)

#ifdef PHIST_HAVE_COLPACK
    call colorcrs(A,2,3,ierr)
#endif

    write(*,*) 'created new crsMat with dimensions', A%nRows, A%nCols, A%nEntries
    flush(6)
    A_ptr = c_loc(A)

    ierr = 0

    !--------------------------------------------------------------------------------
  end subroutine phist_DcrsMat_read_mm


  !==================================================================================
  !> read MatrixMarket file
  subroutine phist_DcrsMat_create_fromRowFunc(A_ptr, comm_ptr,nrows, ncols, maxnne_per_row, rowFunc_ptr, ierr) &
    & bind(C,name='phist_DcrsMat_create_fromRowFunc_f') ! circumvent bug in opari (openmp instrumentalization)
    use, intrinsic :: iso_c_binding
    use env_module, only: newunit
    use mpi
    !--------------------------------------------------------------------------------
    type(C_PTR),        intent(out) :: A_ptr
    type(C_PTR),        value  :: comm_ptr
    integer(kind=C_INT64_T),     value       :: nrows, ncols
    integer(C_INT32_T), value           :: maxnne_per_row
    type(C_FUNPTR),     value       :: rowFunc_ptr
    integer(C_INT),     intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    procedure(matRowFunc), pointer :: rowFunc
    !--------------------------------------------------------------------------------
    integer(kind=G_GIDX_T), allocatable :: idx(:,:)
    real(kind=8), allocatable :: val(:)
    integer(kind=8) :: i, globalRows, globalCols
    integer(kind=G_GIDX_T) :: j, j_, globalEntries
    integer(kind=G_GIDX_T) :: i_
    integer(kind=G_LIDX_T) :: nne
    !integer :: funit
    integer(kind=8) :: localDim(2), globalDim(2)
    real(kind=8) :: wtime
    integer, pointer :: comm
    !--------------------------------------------------------------------------------

    ! get procedure pointer
    call c_f_procpointer(rowFunc_ptr, rowFunc)
    call c_f_pointer(comm_ptr, comm)

    allocate(A)

    ! now read the dimensions
    globalRows = nrows
    globalCols = ncols
    globalEntries = int(maxnne_per_row,kind=8)*int(nrows,kind=8)

    call map_setup(A%row_map, comm, globalRows, ierr)
    if( ierr .ne. 0 ) return
    if( A%row_map%me .eq. 0 ) then
      write(*,*) 'creating matrix from rowFunc'
      write(*,*) 'CrsMat:', globalRows, globalCols, globalEntries
      flush(6)
    end if

    A%nRows = A%row_map%nlocal(A%row_map%me)
    A%nCols = A%nRows
    A%nEntries = int(maxnne_per_row,kind=8)*int(A%nRows,kind=8)

    ! allocate temporary buffers
    allocate(idx(maxnne_per_row,2),&
             val(maxnne_per_row),stat=ierr)

    if (ierr/=0) then
      ierr=-44
      return
    end if

    ! allocate crs matrix
    allocate(A%row_offset(A%nRows+1),&
             A%global_col_idx(A%nEntries),&
             A%val(A%nEntries),stat=ierr)

    if (ierr/=0) then
      ierr=-44
      return
    end if

call mpi_barrier(A%row_map%comm, ierr)
wtime = mpi_wtime()

    ! note: the matFuncs in the physics repo are not
    ! thread-safe (as of Jul 2014), hence the 'ordered' clause.
    ! get data, try to respect NUMA
    A%row_offset(1) = 1_8
!$omp parallel do schedule(static) ordered
    do i = 1, A%nRows, 1
!$omp ordered
      i_ = A%row_map%distrib(A%row_map%me)+i-2
      call rowFunc(i_, nne, idx(:,1), val)
      j = A%row_offset(i)
      j_ = j + int(nne-1,kind=8)
      A%global_col_idx(j:j_) = idx(1:nne,1)+1
      A%val(j:j_) = val(1:nne)
      A%row_offset(i+1) = A%row_offset(i)+int(nne,kind=8)
!$omp end ordered
    end do
    A%nEntries = A%row_offset(A%nRows+1)-1

call mpi_barrier(A%row_map%comm, ierr)
wtime = mpi_wtime() - wtime
if( A%row_map%me .eq. 0 ) then
  write(*,*) 'read matrix from row func in', wtime, 'seconds'
  flush(6)
end if

wtime = mpi_wtime()

    call sort_global_cols(A)

call mpi_barrier(A%row_map%comm, ierr)
wtime = mpi_wtime() - wtime
if( A%row_map%me .eq. 0 ) then
  write(*,*) 'sort global cols in', wtime, 'seconds'
  flush(6)
end if


#ifdef PHIST_HAVE_PARMETIS
wtime = mpi_wtime()

    call repartcrs(A,3,ierr)
    if (ierr/=0) then
      if( A%row_map%me .eq. 0 ) then
        write(*,*) 'repartitioning failed with error code ',ierr
      end if
      return
    end if

call mpi_barrier(A%row_map%comm, ierr)
wtime = mpi_wtime() - wtime
if( A%row_map%me .eq. 0 ) then
  write(*,*) 'repartitioned in', wtime, 'seconds'
end if
#endif
! write matrix to mat.mm
!do i_ = 0, A%row_map%nProcs-1

  !call MPI_Barrier(A%row_map%comm, ierr)
  !if( A%row_map%me .ne. i_ ) cycle

  !if( i_ .eq. 0 ) then
    !open(unit   = newunit(funit), file = 'out.mm', action = 'write')
    !write(funit,'(A)') '%%MatrixMarket matrix coordinate real general'
    !write(funit,'(A)') '% generated by SpinChainSZ'
    !write(funit,*) A%row_map%distrib(A%row_map%nProcs)-1, A%row_map%distrib(A%row_map%nProcs)-1, A%nEntries
  !else
    !open(unit   = newunit(funit), file = 'out.mm', action = 'write', position = 'append')
  !end if

  !do i = 1, A%nRows
    !do j = A%row_offset(i), A%row_offset(i+1)-1, 1
      !write(funit,*) A%row_map%distrib(A%row_map%me)+i-1, A%global_col_idx(j), A%val(j)
    !end do
  !end do
  !close(funit)

!end do

wtime = mpi_wtime()

    call setup_commBuff(A, A%comm_buff)

call mpi_barrier(A%row_map%comm, ierr)
wtime = mpi_wtime() - wtime
if( A%row_map%me .eq. 0 ) then
  write(*,*) 'setup comm buffers in', wtime, 'seconds'
end if


wtime = mpi_wtime()

    call sort_rows_local_nonlocal(A)

call mpi_barrier(A%row_map%comm, ierr)
wtime = mpi_wtime() - wtime
if( A%row_map%me .eq. 0 ) then
  write(*,*) 'sorted rows local/nonlocal in', wtime, 'seconds'
end if

#ifdef PHIST_HAVE_COLPACK
    wtime = mpi_wtime()
    !TODO - make the coloring optional, it only makes sense
    !       for CARP-CG right now
    call colorcrs(A,2,3,ierr)
    if (ierr/=0) then
      if( A%row_map%me .eq. 0 ) then
        write(*,*) 'graph coloring failed'
      end if
      return
    end if
    call mpi_barrier(A%row_map%comm, ierr)
    wtime = mpi_wtime() - wtime
    if( A%row_map%me .eq. 0 ) then
      write(*,*) 'local dist-2 coloring in', wtime, 'seconds'
    end if
#endif

    ! calculate actual global dimension
    localDim(1) = A%nrows
    localDim(2) = A%nEntries
    call mpi_allreduce(localDim, globalDim, 2, MPI_INTEGER8, MPI_SUM, &
      &                A%row_map%comm, ierr)
    if( A%row_map%me .eq. 0 ) then
      write(*,*) 'created new crsMat with dimensions', globalDim
      flush(6)
    end if
    A_ptr = c_loc(A)

    ierr = 0

    !--------------------------------------------------------------------------------
  end subroutine phist_DcrsMat_create_fromRowFunc


  !==================================================================================
  subroutine phist_DcrsMat_delete(A_ptr, ierr) bind(C,name='phist_DcrsMat_delete_f')
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),    value       :: A_ptr
    integer(C_INT), intent(out) :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    integer :: i
#if defined(TESTING) && PHIST_OUTLEV >= 4
    integer(C_INTPTR_T) :: dummy
#endif
    !--------------------------------------------------------------------------------

#if defined(TESTING) && PHIST_OUTLEV >= 4
    write(*,*) 'deleting crsMat at address', transfer(A_ptr,dummy)
    flush(6)
#endif
    if( c_associated(A_ptr) ) then
      ! get pointer
      call c_f_pointer(A_ptr,A)

      ! free all mpi requests
      if( allocated(A%comm_buff%recvData) ) then
        call delete_buffers(A%comm_buff)
      end if

      ! deallocate matrix
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
    
    ierr=0

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

    call crsmat_times_mvec(alpha,A,shifts,x,beta,y,ierr)

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

    call crsmat_times_mvec(alpha,A,shifts,x,beta,y,ierr)

  end subroutine phist_DcrsMat_times_mvec_vadd_mvec

  !==================================================================================

  subroutine phist_Dcarp_setup(A_ptr, numShifts, shifts_r, shifts_i, nrms_ai2i_ptr, work_ptr,ierr) &
    & bind(C,name='phist_Dcarp_setup_f') ! circumvent bug in opari (openmp instrumentalization)
    use, intrinsic :: iso_c_binding
    !--------------------------------------------------------------------------------
    type(C_PTR),      value         :: A_ptr
    integer(C_INT),   value         :: numShifts
    real(kind=c_double), intent(in) :: shifts_r(numShifts), shifts_i(numShifts)
    TYPE(C_PTR), value :: nrms_ai2i_ptr
    real(kind=8), pointer, dimension(:,:) :: nrms_ai2i
    type(C_PTR), intent(out) :: work_ptr
    integer(C_INT),   intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    ! this array contains the 1/np, where np is the number of procs contributing
    ! to a value in the 'receive' area (the elements sent during regular spMVM)
    type(sparseArray_t), pointer :: invProcCount

    !--------------------------------------------------------------------------------
    
    integer, allocatable, dimension(:) :: procCount !! temporary array
    integer :: theShape(2)
    integer :: i,k,l
    integer(kind=8) :: j

    if( .not. c_associated(A_ptr)) then
      ierr = -88
      return
    end if
!write(*,*) 'numShifts=',numShifts
!do i=1,numShifts
!  write(*,*) shifts_r(i),' + i',shifts_i(i)
!end do

    call c_f_pointer(A_ptr,A)
    theShape(1)=A%nRows
    theShape(2)=numShifts
    call c_f_pointer(nrms_ai2i_ptr,nrms_ai2i,theShape)

    ! start by putting the diagonal elements of A in the first column
    ! of this array, will be overwritten by the kernel

!$omp parallel do private(j) schedule(static)
    do i = 1,A%nRows
      nrms_ai2i(i,:)=0.d0
      do j = A%row_offset(i), A%nonlocal_offset(i)-1, 1
        if (A%col_idx(j).eq.i) then
          nrms_ai2i(i,1)=A%val(j)
        end if
      end do
    end do

    ! compute inverse row norms for shift i in column i
    call crsmat_norms_ai2i(numShifts, A%nRows, A%nEntries, &
        A%row_offset, A%val, shifts_r,shifts_i,nrms_ai2i)
    
    ! work is used to store the averaging coefficients for the x halo
    ! received from other processes (in a sparseArray_t object)
    allocate(procCount(A%nRows),stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    procCount(1:A%nRows)=1 ! contains for every owned node the number
                           ! of contributing procs (default 1, just me)
    
    do i=1,A%comm_buff%nSendProcs, 1
      k = A%comm_buff%sendInd(i)
      l = A%comm_buff%sendInd(i+1)
      do j = k, l-1, 1
        procCount(A%comm_buff%sendRowBlkInd(j)) = &
        procCount(A%comm_buff%sendRowBlkInd(j)) + 1
      end do
    end do

    allocate(invProcCount)
    invProcCount%n=A%nRows
    invProcCount%nnz=count(procCount.gt.1)
    allocate(invProcCount%idx(invProcCount%nnz),&
             invProcCount%val(invProcCount%nnz), stat=ierr)
    if (ierr/=0) then
      ierr=-44
      return
    end if
    j=0
    do i=1,A%nRows
      if (procCount(i) .gt. 1) then
        j=j+1
        invProcCount%val(j)=1.0/dble(procCount(i))
        invProcCount%idx(j)=i
      end if
    end do
    deallocate(procCount)
    work_ptr=c_loc(invProcCount)    
    
    !if (A%row_map%me==0 .and. A%row_map%coloringType==2) then
      !write(6,*) 'CARP will exploit local dist-2 coloring'
      !flush(6)
    !end if
  end subroutine phist_Dcarp_setup

  subroutine phist_Dcarp_sweep(A_ptr, numShifts, shifts_r, shifts_i, &
        b_ptr, x_r_ptr, x_i_ptr, nrms_ai2i_ptr, work_ptr, omegas, ierr) &
  bind(C,name='phist_Dcarp_sweep_f')
    use, intrinsic :: iso_c_binding
    implicit none
    !--------------------------------------------------------------------------------
    type(C_PTR),      value         :: A_ptr, b_ptr
    integer(c_int),   value    :: numShifts
    type(C_PTR)                     :: x_r_ptr(numShifts), x_i_ptr(numShifts)
    real(kind=c_double), intent(in) :: shifts_r(numShifts), shifts_i(numShifts)
    real(kind=c_double), intent(in) :: omegas(numShifts)
    type(C_PTR),      value         :: nrms_ai2i_ptr,work_ptr
    integer(C_INT),   intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    real(kind=8), pointer :: nrms_ai2i(:,:)
    type(MVec_t), pointer :: x_r, x_i, b
    type(sparseArray_t), pointer :: invProcCount
    
    !--------------------------------------------------------------------------------
    integer :: iSys
    integer :: ldx, ldb, nvec
    integer :: theShape(2)
    integer :: sendBuffSize,recvBuffSize
    logical :: b_is_zero
    real(kind=8), dimension(0,0), target :: bzero

    ierr=0
    
    if ( .not. c_associated(A_ptr) ) then
      write(*,*) 'A is NULL'
      ierr = -88
      return
    else
      call c_f_pointer(A_ptr,A)
    end if

    if ( .not. c_associated(b_ptr) ) then
      b_is_zero=.true.
      ldb=0
    else
      b_is_zero=.false.
      call c_f_pointer(b_ptr,b)
      ldb = size(b%val,1)
    end if

    theShape(1) = A%nRows
    theShape(2) = numShifts

    !write(*,*) 'enter carp_sweep_f, bzero=',b_is_zero
    if ( .not. c_associated(nrms_ai2i_ptr) ) then
      write(*,*) 'nrms_ai2i is NULL'
      ierr = -88
      return
    else
      call c_f_pointer(nrms_ai2i_ptr,nrms_ai2i,theShape)
    end if

    if ( .not. c_associated(work_ptr) ) then
      write(*,*) 'work is NULL'
      ierr = -88
      return
    else
      call c_f_pointer(work_ptr,invProcCount)
    end if

    ! treat one shift at a time for the moment, here there
    ! is potential for additional parallelism, of course,
    ! but the user can also handle this level himself by
    ! passsing in one shift at a time to this function.
    do iSys=1,numShifts
      
      ! check that all C pointers are non-null
      if ( .not. c_associated(x_r_ptr(iSys)) .or. &
         & .not. c_associated(x_i_ptr(iSys)) ) then
        write(*,*) 'an input mvec x is NULL, index ',iSys
        ierr = -88
        return
      end if
      call c_f_pointer(x_r_ptr(iSys),x_r)
      call c_f_pointer(x_i_ptr(iSys),x_i)

      nvec = x_r%jmax-x_r%jmin+1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! (re-)allocate communication buffers,
      ! space for 2*nvecs vector halos is needed
      ! because we have x_r and x_i (complex x).
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      if( allocated(A%comm_buff%recvData) ) then
        if( size(A%comm_buff%recvData,1) .ne. nvec ) then
          call delete_buffers(A%comm_buff)
        end if
      end if

      sendBuffSize = A%comm_buff%sendInd(A%comm_buff%nSendProcs+1)-1
      recvBuffSize = A%comm_buff%recvInd(A%comm_buff%nRecvProcs+1)-1

      if( .not. allocated(A%comm_buff%recvData) ) then
        call alloc_buffers(A%comm_buff, A%row_map%comm, nvec,ierr)
        if (ierr/=0) then
          ierr=-44
          return
        end if
      end if

      ldx = size(x_r%val,1)
    
      !write(*,*) 'nvec=',nvec
      !write(*,*) 'x_r%jmin=',x_r%jmin
      !write(*,*) 'x_r%jmax=',x_r%jmax
      !write(*,*) 'bounds(1)=',lbound(x_r%val,1),ubound(x_r%val,1)
      !write(*,*) 'bounds(2)=',lbound(x_r%val,2),ubound(x_r%val,2)
      !write(*,*) 'ldx=',ldx
    
      ! exchange necessary elements of X. This is the same      
      ! operation as for an spMVM: import from the domain       
      ! map into the column map of A.                           
      call Dcarp_import2(A,x_r,x_i,ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_import returned ierr=',ierr
        return
      end if

      ! apply Kaczmarz forward sweep
      if (b_is_zero) then
      call dkacz_selector(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys),shifts_i(iSys), &
                bzero, ldb, &
                x_r%val(x_r%jmin,1),x_i%val(x_i%jmin,1), ldx, &
                A%comm_buff%recvData(1,1,1),&
                A%comm_buff%recvData(1,1,2),&
                nrms_ai2i(1,iSys),omegas(iSys),&
                1,A%nRows,+1)
        else
      call dkacz_selector(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys),shifts_i(iSys), &
                b%val(b%jmin,1), ldb, &
                x_r%val(x_r%jmin,1),x_i%val(x_i%jmin,1), ldx, &
                A%comm_buff%recvData(1,1,1),&
                A%comm_buff%recvData(1,1,2),&
                nrms_ai2i(1,iSys),omegas(iSys),&
                1,A%nRows,+1)        
        end if
      ! exchange values and average (export/average operation
      ! from col map into domain map)
      call Dcarp_average2(A,x_r, x_i, invProcCount, ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_average returned ierr=',ierr
        return
      end if
      ! import again into the column map (get updated halo elements)
      ! before the backward sweep
      call Dcarp_import2(A,x_r,x_i,ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_import returned ierr=',ierr
        return
      end if

      ! apply Kaczmarz backward sweep
      !write(*,*) 'kacz (b)'
      if (b_is_zero) then
      call dkacz_selector(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys),shifts_i(iSys), &
                bzero, ldb, &
                x_r%val(x_r%jmin,1),x_i%val(x_i%jmin,1), ldx, &
                A%comm_buff%recvData(1,1,1),&
                A%comm_buff%recvData(1,1,2),&
                nrms_ai2i(1,iSys),omegas(iSys),&
                A%nRows,1,-1)      
      else
      call dkacz_selector(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys),shifts_i(iSys), &
                b%val(b%jmin,1), ldb, &
                x_r%val(x_r%jmin,1),x_i%val(x_i%jmin,1), ldx, &
                A%comm_buff%recvData(1,1,1),&
                A%comm_buff%recvData(1,1,2),&
                nrms_ai2i(1,iSys),omegas(iSys),&
                A%nRows,1,-1)
      end if
      ! exchange values and average, stay in domain map
      ! (i.e. do not import halo until the next call to
      ! carp_sweep)
      call Dcarp_average2(A,x_r, x_i, invProcCount, ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_average returned ierr=',ierr
        return
      end if
    end do ! shifts
    
  end subroutine phist_Dcarp_sweep


  ! variant of carp_sweep with real matrix, real shift, and real vectors
  subroutine phist_Dcarp_sweep_real(A_ptr, numShifts, shifts_r, &
        b_ptr, x_r_ptr, nrms_ai2i_ptr, work_ptr, omegas, ierr) &
  bind(C,name='phist_Dcarp_sweep_real_f')
    use, intrinsic :: iso_c_binding
    implicit none
    !--------------------------------------------------------------------------------
    type(C_PTR),      value         :: A_ptr, b_ptr
    integer(c_int),   value    :: numShifts
    type(C_PTR)                     :: x_r_ptr(numShifts)
    real(kind=c_double), intent(in) :: shifts_r(numShifts)
    real(kind=c_double), intent(in) :: omegas(numShifts)
    type(C_PTR),      value         :: nrms_ai2i_ptr,work_ptr
    integer(C_INT),   intent(out)   :: ierr
    !--------------------------------------------------------------------------------
    type(CrsMat_t), pointer :: A
    real(kind=8), pointer :: nrms_ai2i(:,:)
    type(MVec_t), pointer :: x_r, b
    type(sparseArray_t), pointer :: invProcCount
    
    !--------------------------------------------------------------------------------
    integer :: iSys
    integer :: ldx, ldb, nvec
    integer :: theShape(2)
    integer :: sendBuffSize,recvBuffSize
    logical :: b_is_zero
    real(kind=8), dimension(0,0), target :: bzero

    ierr=0
    
    if ( .not. c_associated(A_ptr) ) then
      write(*,*) 'A is NULL'
      ierr = -88
      return
    else
      call c_f_pointer(A_ptr,A)
    end if

    if ( .not. c_associated(b_ptr) ) then
      b_is_zero=.true.
      ldb=0
    else
      b_is_zero=.false.
      call c_f_pointer(b_ptr,b)
      ldb = size(b%val,1)
    end if

    theShape(1) = A%nRows
    theShape(2) = numShifts

    !write(*,*) 'enter carp_sweep_f, bzero=',b_is_zero
    if ( .not. c_associated(nrms_ai2i_ptr) ) then
      write(*,*) 'nrms_ai2i is NULL'
      ierr = -88
      return
    else
      call c_f_pointer(nrms_ai2i_ptr,nrms_ai2i,theShape)
    end if

    if ( .not. c_associated(work_ptr) ) then
      write(*,*) 'work is NULL'
      ierr = -88
      return
    else
      call c_f_pointer(work_ptr,invProcCount)
    end if

    ! treat one shift at a time for the moment, here there
    ! is potential for additional parallelism, of course,
    ! but the user can also handle this level himself by
    ! passsing in one shift at a time to this function.
    do iSys=1,numShifts
      
      ! check that all C pointers are non-null
      if ( .not. c_associated(x_r_ptr(iSys)) ) then
        write(*,*) 'an input mvec x is NULL, index ',iSys
        ierr = -88
        return
      end if
      call c_f_pointer(x_r_ptr(iSys),x_r)

      nvec = x_r%jmax-x_r%jmin+1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! (re-)allocate communication buffers,
      ! space for 2*nvecs vector halos is needed
      ! because we have x_r and x_i (complex x).
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      if( allocated(A%comm_buff%recvData) ) then
        if( size(A%comm_buff%recvData,1) .ne. nvec ) then
          call delete_buffers(A%comm_buff)
        end if
      end if

      sendBuffSize = A%comm_buff%sendInd(A%comm_buff%nSendProcs+1)-1
      recvBuffSize = A%comm_buff%recvInd(A%comm_buff%nRecvProcs+1)-1

      if( .not. allocated(A%comm_buff%recvData) ) then
        call alloc_buffers(A%comm_buff, A%row_map%comm, nvec,ierr)
        if (ierr/=0) then
          ierr=-44
          return
        end if
      end if

      ldx = size(x_r%val,1)
    
      !write(*,*) 'nvec=',nvec
      !write(*,*) 'x_r%jmin=',x_r%jmin
      !write(*,*) 'x_r%jmax=',x_r%jmax
      !write(*,*) 'bounds(1)=',lbound(x_r%val,1),ubound(x_r%val,1)
      !write(*,*) 'bounds(2)=',lbound(x_r%val,2),ubound(x_r%val,2)
      !write(*,*) 'ldx=',ldx
    
      ! exchange necessary elements of X. This is the same      
      ! operation as for an spMVM: import from the domain       
      ! map into the column map of A.                           
      call Dcarp_import(A,x_r,ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_import returned ierr=',ierr
        return
      end if

      ! apply Kaczmarz forward sweep
      if (b_is_zero) then
      call dkacz_selector_real(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys), &
                bzero, ldb, &
                x_r%val(x_r%jmin,1),ldx, &
                A%comm_buff%recvData,&
                nrms_ai2i(1,iSys),omegas(iSys),&
                1,A%nRows,+1)
        else
      call dkacz_selector_real(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys), &
                b%val(b%jmin,1), ldb, &
                x_r%val(x_r%jmin,1), ldx, &
                A%comm_buff%recvData,&
                nrms_ai2i(1,iSys),omegas(iSys),&
                1,A%nRows,+1)        
        end if
      ! exchange values and average (export/average operation
      ! from col map into domain map)
      call Dcarp_average(A,x_r, invProcCount, ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_average returned ierr=',ierr
        return
      end if
      ! import again into the column map (get updated halo elements)
      ! before the backward sweep
      call Dcarp_import(A,x_r,ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_import returned ierr=',ierr
        return
      end if

      ! apply Kaczmarz backward sweep
      !write(*,*) 'kacz (b)'
      if (b_is_zero) then
      call dkacz_selector_real(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys), &
                bzero, ldb, &
                x_r%val(x_r%jmin,1),ldx, &
                A%comm_buff%recvData,&
                nrms_ai2i(1,iSys),omegas(iSys),&
                A%nRows,1,-1)      
      else
      call dkacz_selector_real(nvec, A%nRows, recvBuffSize,A%nCols, A%nEntries, &
                A%row_offset, A%nonlocal_offset, A%col_idx, A%val, &
                A%row_map,&
                shifts_r(iSys), &
                b%val(b%jmin,1), ldb, &
                x_r%val(x_r%jmin,1), ldx, &
                A%comm_buff%recvData,&
                nrms_ai2i(1,iSys),omegas(iSys),&
                A%nRows,1,-1)
      end if
      ! exchange values and average, stay in domain map
      ! (i.e. do not import halo until the next call to
      ! carp_sweep)
      call Dcarp_average(A,x_r, invProcCount, ierr)
      if (ierr/=0) then
        write(*,*) 'call Dcarp_average returned ierr=',ierr
        return
      end if
    end do ! shifts
    
  end subroutine phist_Dcarp_sweep_real

  subroutine phist_Dcarp_destroy(A_ptr, numShifts, work_ptr, ierr) &
  bind(C,name='phist_Dcarp_destroy_f')
    use, intrinsic :: iso_c_binding

    type(C_PTR), value :: A_ptr
    integer(c_int), value :: numShifts
    type(C_PTR), value :: work_ptr
    integer(c_int), intent(out) :: ierr
    !-----------------------------------------!
    TYPE(sparseArray_t), pointer :: invProcCount
    
    ! the C top layer allocated nrms_ai2i, so it also deletes them
    ierr=0
    
    call c_f_pointer(work_ptr, invProcCount)
    if (allocated(invProcCount%idx)) then
      deallocate(invProcCount%idx)
    end if
    if (allocated(invProcCount%val)) then
      deallocate(invProcCount%val)
    end if
    
  end subroutine phist_Dcarp_destroy

subroutine Dspmv_import(A,x,ierr)
use mpi
implicit none

  TYPE(crsMat_t) :: A
  TYPE(mvec_t)   :: x
  integer, intent(out) :: ierr
  ! local
  integer :: nvec, ldx, sendBuffSize
  logical :: strided_x
  integer :: i,j,k,l

    nvec = x%jmax-x%jmin+1
    ldx = size(x%val,1)

    ! determin data layout
    if( .not. x%is_view .or. &
      & ( x%jmin .eq. lbound(x%val,1) .and. &
      &   x%jmax .eq. ubound(x%val,1)       ) ) then
      strided_x = .false.
    else
      strided_x = .true.
    end if


    ! exchange necessary elements
!write(*,*) 'CRS', A%row_map%me, 'sendRowBlkInd', A%comm_buff%sendRowBlkInd
!write(*,*) 'CRS', A%row_map%me, 'recvRowBlkInd', A%comm_buff%recvRowBlkInd

    ! start buffer irecv
    ! we could also set up persistent communication channels here... and use MPI_Startall later
    if( allocated(A%comm_buff%recvData) ) then
      if( size(A%comm_buff%recvData,1) .ne. nvec ) then
        call delete_buffers(A%comm_buff)
      end if
    end if
    if( .not. allocated(A%comm_buff%recvData) ) then
      call alloc_buffers(A%comm_buff,A%row_map%comm,nvec,ierr)
    end if

    ! start receiving (e.g. irecv)
    if( A%comm_buff%nRecvProcs .gt. 0 ) then
      !call mpi_startall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,A%comm_buff%recvStatus,ierr)
      call mpi_startall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,ierr)
    end if


    sendBuffSize = size(A%comm_buff%sendBuffInd,1)
    ! gather send data in buffer
    if( nvec .eq. 1 ) then
      if( strided_x ) then
        call dspmv_buff_cpy_strided_1(A%nrows, sendBuffSize, ldx, x%val(x%jmin,1), A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      else
        call dspmv_buff_cpy_1(A%nrows, sendBuffSize, x%val, A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      end if
    else if( nvec .eq. 2 ) then
      if( strided_x ) then
        call dspmv_buff_cpy_strided_2(A%nrows, sendBuffSize, ldx, x%val(x%jmin,1), A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      else
        call dspmv_buff_cpy_2(A%nrows, sendBuffSize, x%val, A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      end if
    else if( nvec .eq. 4 ) then
      if( strided_x ) then
        call dspmv_buff_cpy_strided_4(A%nrows, sendBuffSize, ldx, x%val(x%jmin,1), A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      else
        call dspmv_buff_cpy_4(A%nrows, sendBuffSize, x%val, A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      end if
    else if( nvec .eq. 8 ) then
      if( strided_x ) then
        call dspmv_buff_cpy_strided_8(A%nrows, sendBuffSize, ldx, x%val(x%jmin,1), A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      else
        call dspmv_buff_cpy_8(A%nrows, sendBuffSize, x%val, A%comm_buff%sendBuffInd, A%comm_buff%sendData)
      end if
    else
      call dspmv_buff_cpy_general(nvec, A%nrows, sendBuffSize, ldx, x%val(x%jmin,1), A%comm_buff%sendBuffInd, A%comm_buff%sendData)
    end if


    ! start sending (e.g. issend)
    if( A%comm_buff%nSendProcs .gt. 0 ) then
      !call mpi_startall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
      call mpi_startall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,ierr)
    end if

    ! wait till all data arrived
    if( A%comm_buff%nRecvProcs .gt. 0 ) then
      call mpi_waitall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,A%comm_buff%recvStatus,ierr)
    end if

    ! also wait till all data was sent to circumvent dumb MPI implementations that don't do that in the background!
    if( A%comm_buff%nSendProcs .gt. 0 ) then
      call mpi_waitall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
    end if

  end subroutine Dspmv_import

! import halo elements before a CARP sweep. This is the same
! operation performed before an spMVM. On exit, the halo elements
! are found in A%comm_buf%recvData(1:nvec,:,1)
subroutine Dcarp_import(A,x, ierr)
  use mpi
  implicit none

  TYPE(crsMat_t) :: A
  TYPE(mvec_t)   :: x
  integer, intent(out) :: ierr
  ! local
  integer :: nvec,nrecv
  integer :: i,j,k,l

  nvec=x%jmax-x%jmin+1
  nrecv=size(A%comm_buff%recvData,1)
  if (nrecv/=nvec) then
    write(*,*) 'WARNING: unexpected buffer size in ',__FILE__,' line ',__LINE__
  end if

      
  ! start buffer irecv
    if( A%comm_buff%nRecvProcs .gt. 0 ) then
      call mpi_startall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,ierr)
    end if

  ! start buffer isend
  do i=1,A%comm_buff%nSendProcs, 1
    k = A%comm_buff%sendInd(i)
    l = A%comm_buff%sendInd(i+1)
    do j = k, l-1, 1
      A%comm_buff%sendData(1:nvec,j,1) = &
        x%val(x%jmin:x%jmax,A%comm_buff%sendRowBlkInd(j))
    end do
  end do

    ! start sending (e.g. issend)
    if( A%comm_buff%nSendProcs .gt. 0 ) then
      !call mpi_startall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
      call mpi_startall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,ierr)
    end if

  if( A%comm_buff%nRecvProcs .gt. 0 ) then
    ! wait till all data arrived
    call mpi_waitall(A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,A%comm_buff%recvStatus,ierr)
  end if

  if( A%comm_buff%nSendProcs .gt. 0 ) then
    ! wait till all data was sent
    call mpi_waitall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
  end if
end subroutine Dcarp_import

! like carp_import, but for two vector blocks. Halo elements are put into
! A%comm_buff%recvData(1:nvec,:,K) for xK, K=1,2
subroutine Dcarp_import2(A,x1,x2, ierr)
  use mpi
  implicit none

  TYPE(crsMat_t) :: A
  TYPE(mvec_t)   :: x1,x2
  integer, intent(out) :: ierr

  ! local
  integer :: nvec,nrecv
  integer :: i,j,k,l

  nvec=x1%jmax-x1%jmin+1
  if (nvec/=x2%jmax-x2%jmin+1) then
    ierr = -55
    return
  end if
  nrecv=size(A%comm_buff%recvData,1)
  if (nrecv/=nvec) then
    ! this should have been adjusted (if necessary) beforehand
    ierr=-55
    return
  end if

      
  ! start buffer irecv
    if( A%comm_buff%nRecvProcs .gt. 0 ) then
      call mpi_startall(2*A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,ierr)
    end if

  ! start buffer isend
  do i=1,A%comm_buff%nSendProcs, 1
    k = A%comm_buff%sendInd(i)
    l = A%comm_buff%sendInd(i+1)
    do j = k, l-1, 1
      A%comm_buff%sendData(1:nvec,j,1) = &
        x1%val(x1%jmin:x1%jmax,A%comm_buff%sendRowBlkInd(j))
      A%comm_buff%sendData(1:nvec,j,2) = &
        x2%val(x2%jmin:x2%jmax,A%comm_buff%sendRowBlkInd(j))
    end do
  end do

    ! start sending (e.g. issend)
    if( A%comm_buff%nSendProcs .gt. 0 ) then
      !call mpi_startall(A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
      call mpi_startall(2*A%comm_buff%nSendProcs,A%comm_buff%sendRequests,ierr)
    end if

  if( A%comm_buff%nRecvProcs .gt. 0 ) then
    ! wait till all data arrived
    call mpi_waitall(2*A%comm_buff%nRecvProcs,A%comm_buff%recvRequests,A%comm_buff%recvStatus,ierr)
  end if

  if( A%comm_buff%nSendProcs .gt. 0 ) then
    ! wait till all data was sent
    call mpi_waitall(2*A%comm_buff%nSendProcs,A%comm_buff%sendRequests,A%comm_buff%sendStatus,ierr)
  end if
end subroutine Dcarp_import2

! private subroutine to exchange halo elements and average
!
! The standard spMVM import is: get values of halo elements
! from neighbors (import vector from domain map to column  
! map). We now do the opposite: we receive alternative values
! for the local nodes ghosted on other procs and update them
! by averaging. Before the next spMVM or CARP sweep, a standard
! halo import occurs to update the values on the other procs.
subroutine Dcarp_average(A,x,invProcCount, ierr)
  use mpi
  implicit none

  TYPE(crsMat_t) :: A
  TYPE(mvec_t)   :: x
  ! this array keeps track from how many procs we receive
  ! contributions to a local x entry. In principle we could
  ! set this up once and for all.
  ! (TODO, or maybe it already exists somewhere?)
  type(sparseArray_t), intent(in)  :: invProcCount
  integer :: ierr

  ! local
  integer :: nvec
  integer :: i,j,k,l

  nvec=x%jmax-x%jmin+1

  ! in this subroutine we assume that the (updated) halo elements
  ! are still in the receive buffer (recvData) from the previous 
  ! import in carp_sweep.

    ! NOTE: in the following communication step, send and recv bufs are deliberately 
    !       swapped! This is because we now *send* our updated halo elements and *recv*
    !       alternative values for local nodes, which we average.

      ! start buffered *irecv* of clones of local values
      if (A%comm_buff%nSendProcs .gt. 0) then
        call mpi_startall(A%comm_buff%nSendProcs,A%comm_buff%inv_recvRequests,ierr)
      end if

      ! start buffered *isend* of my results in the halo
      ! NOTE: we assume that the data we're sending is already
      ! in the recvBuf. This is the case in carp_sweep: the halo
      ! is imported and the kernel which is called (in kacz_selector)
      ! operates directly on the recvBuf and updates its contents.
      if (A%comm_buff%nRecvProcs .gt. 0) then
        call mpi_startall(A%comm_buff%nRecvProcs,A%comm_buff%inv_sendRequests,ierr)
      end if

    if( A%comm_buff%nSendProcs .gt. 0 ) then
      ! wait till all data arrived
      call mpi_waitall(A%comm_buff%nSendProcs,&
                       A%comm_buff%inv_recvRequests,&
                       A%comm_buff%sendStatus,ierr)
    end if
    
    ! now average the values in the local elements for which we have received
    ! alternative values from the neighbors
    do i=1,A%comm_buff%nSendProcs
      k = A%comm_buff%sendInd(i)
      l = A%comm_buff%sendInd(i+1)
      do j = k, l-1, 1

          x%val(x%jmin:x%jmax,A%comm_buff%sendRowBlkInd(j)) &
        = x%val(x%jmin:x%jmax,A%comm_buff%sendRowBlkInd(j)) &
        + A%comm_buff%sendData(1:nvec,j,1)
      end do
    end do

    ! average: divide entries summed up from others by the number of
    ! contributors (in invProcCount constructed in carp_setup)
    do i=1,invProcCount%nnz
      k=invProcCount%idx(i)
      x%val(x%jmin:x%jmax,k) = x%val(x%jmin:x%jmax,k)*invProcCount%val(i)
    end do

    if( A%comm_buff%nSendProcs .gt. 0 ) then
      ! wait till all data was sent
      call mpi_waitall(A%comm_buff%nRecvProcs,&
                       A%comm_buff%inv_sendRequests,&
                       A%comm_buff%recvStatus,ierr)
    end if

end subroutine Dcarp_average

subroutine Dcarp_average2(A,x1,x2,invProcCount, ierr)
  use mpi
  implicit none

  TYPE(crsMat_t) :: A
  TYPE(mvec_t)   :: x1,x2
  ! this array keeps track from how many procs we receive
  ! contributions to a local x entry. In principle we could
  ! set this up once and for all.
  ! (TODO, or maybe it already exists somewhere?)
  type(sparseArray_t), intent(in)  :: invProcCount
  integer :: ierr

  ! local
  integer :: nvec
  integer :: i,j,k,l

  nvec=x1%jmax-x1%jmin+1
  if (nvec/=x2%jmax-x2%jmin+1) then
    ierr=-55
    return
  end if

  ! in this subroutine we assume that the (updated) halo elements
  ! are still in the receive buffer (recvData) from the previous 
  ! import in carp_sweep.

    ! NOTE: in the following communication step, send and recv bufs are deliberately 
    !       swapped! This is because we now *send* our updated halo elements and *recv*
    !       alternative values for local nodes, which we average.

      ! start buffered *irecv* of clones of local values
      if (A%comm_buff%nSendProcs .gt. 0) then
        call mpi_startall(2*A%comm_buff%nSendProcs,A%comm_buff%inv_recvRequests,ierr)
      end if

      ! start buffered *isend* of my results in the halo
      ! NOTE: we assume that the data we're sending is already
      ! in the recvBuf. This is the case in carp_sweep: the halo
      ! is imported and the kernel which is called (in kacz_selector)
      ! operates directly on the recvBuf and updates its contents.
      if (A%comm_buff%nRecvProcs .gt. 0) then
        call mpi_startall(2*A%comm_buff%nRecvProcs,A%comm_buff%inv_sendRequests,ierr)
      end if

    if( A%comm_buff%nSendProcs .gt. 0 ) then
      ! wait till all data arrived
      call mpi_waitall(2*A%comm_buff%nSendProcs,&
                       A%comm_buff%inv_recvRequests,&
                       A%comm_buff%sendStatus,ierr)
    end if
    
    ! now average the values in the local elements for which we have received
    ! alternative values from the neighbors
    do i=1,A%comm_buff%nSendProcs
      k = A%comm_buff%sendInd(i)
      l = A%comm_buff%sendInd(i+1)
      do j = k, l-1, 1

          x1%val(x1%jmin:x1%jmax,A%comm_buff%sendRowBlkInd(j)) &
        = x1%val(x1%jmin:x1%jmax,A%comm_buff%sendRowBlkInd(j)) &
        + A%comm_buff%sendData(1:nvec,j,1)

          x2%val(x2%jmin:x2%jmax,A%comm_buff%sendRowBlkInd(j)) &
        = x2%val(x2%jmin:x2%jmax,A%comm_buff%sendRowBlkInd(j)) &
        + A%comm_buff%sendData(1:nvec,j,2)
      end do
    end do

    ! average: divide entries summed up from others by the number of
    ! contributors (in invProcCount constructed in carp_setup)
    do i=1,invProcCount%nnz
      k=invProcCount%idx(i)
      x1%val(x1%jmin:x1%jmax,k) = x1%val(x1%jmin:x1%jmax,k)*invProcCount%val(i)
      x2%val(x1%jmin:x2%jmax,k) = x2%val(x2%jmin:x2%jmax,k)*invProcCount%val(i)
    end do

    if( A%comm_buff%nSendProcs .gt. 0 ) then
      ! wait till all data was sent
      call mpi_waitall(2*A%comm_buff%nRecvProcs,&
                       A%comm_buff%inv_sendRequests,&
                       A%comm_buff%recvStatus,ierr)
    end if


end subroutine Dcarp_average2

end module crsmat_module

