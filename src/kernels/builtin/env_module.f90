!> \file env_module.f90
!! Defines env_module, auxiliary functions for the phist builtin kernels
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

!> Implements phist_comm_* subroutines for the builtin kernels and some auxiliary functions
module env_module
  implicit none
  private


  public :: newunit
  public :: phist_comm_create
  public :: phist_comm_delete
  public :: phist_comm_get_rank
  public :: phist_comm_get_size
  !public :: init_random_seed

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


  !================================================================================
  ! seed fortran "random_number", copied from gcc.gnu.org
  subroutine init_random_seed() bind(C,name="init_random_seed")
    use iso_fortran_env, only: int64
#ifdef __INTEL_COMPILER
    use ifport, only: getpid
#endif
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
  
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed




end module env_module

! just used to prevent compiler warnings about unused variables
subroutine TOUCH_EMPTY_FUNCTION(x)
end subroutine TOUCH_EMPTY_FUNCTION
