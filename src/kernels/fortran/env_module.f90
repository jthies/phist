module env_module
  implicit none
  private


  public :: newunit
  public :: phist_comm_create
  public :: phist_comm_delete
  public :: phist_comm_get_rank
  public :: phist_comm_get_size

contains

  !================================================================================
  !> returns an unused unit for open
  !! use with open(unit=newunit(myunit), ...)
  !! (from fortranwiki.org)
  function newunit(aunit)
    !------------------------------------------------------------
    integer                        :: newunit
    integer, intent(out), optional :: aunit
    !------------------------------------------------------------
    integer, parameter :: LUN_MIN = 200, LUN_MAX = 1000
    logical            :: isopen
    integer            :: lun
    !------------------------------------------------------------

    ! search for unused unit
    newunit = -1
    do lun = LUN_MIN, LUN_MAX
      inquire(unit=lun,opened=isopen)
      if( .not. isopen ) then
        newunit = lun
        exit
      end if
    end do

    if( present(aunit) ) aunit = newunit

  end function newunit


  !================================================================================
  ! some simple wrapper routines

  subroutine phist_comm_create(comm_ptr, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    intent(out) :: comm_ptr
    integer(C_INT), intent(out) :: ierr
    !------------------------------------------------------------
    integer, pointer :: comm
    !------------------------------------------------------------
    allocate(comm)
    comm = MPI_COMM_WORLD
    comm_ptr = c_loc(comm)
    ierr = 0
  end subroutine phist_comm_create


  subroutine phist_comm_delete(comm_ptr,ierr) bind(C)
    use, intrinsic :: iso_c_binding
    !------------------------------------------------------------
    type(C_PTR),    value       :: comm_ptr
    integer(C_iNT), intent(out) :: ierr
    !------------------------------------------------------------
    integer, pointer :: comm
    !------------------------------------------------------------
    call c_f_pointer(comm_ptr,comm)
    deallocate(comm)
    ierr = 0
  end subroutine phist_comm_delete

  subroutine phist_comm_get_rank(comm_ptr, rank, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    value       :: comm_ptr
    integer(C_INT), intent(out) :: rank, ierr
    !------------------------------------------------------------
    integer, pointer :: comm
    !------------------------------------------------------------
    call c_f_pointer(comm_ptr,comm)
    call mpi_comm_rank(comm,rank,ierr)
  end subroutine phist_comm_get_rank

  subroutine phist_comm_get_size(comm_ptr, nprocs, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    value       :: comm_ptr
    integer(C_INT), intent(out) :: nprocs, ierr
    !------------------------------------------------------------
    integer, pointer :: comm
    !------------------------------------------------------------
    call c_f_pointer(comm_ptr,comm)
    call mpi_comm_size(comm,nprocs,ierr)
  end subroutine phist_comm_get_size


end module env_module
