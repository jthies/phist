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
  public :: allocate_aligned
  public :: deallocate_aligned
  !public :: init_random_seed

  !================================================================================
  !> helper function to allocate an array of aligned memory
  interface allocate_aligned
    module procedure allocate_aligned_1
    module procedure allocate_aligned_2
    module procedure allocate_aligned_3
  end interface allocate_aligned

  ! C declarations
  interface
    function posix_memalign(p, align, n) bind(C)
      use, intrinsic :: iso_c_binding, only: C_PTR, C_SIZE_T, C_INT
      integer(kind=C_SIZE_T), value :: align, n
      type(C_PTR), intent(out) :: p
      integer(kind=C_INT) :: posix_memalign
    end function posix_memalign
    subroutine free(p) bind(C)
      use, intrinsic :: iso_c_binding, only: C_PTR
      type(C_PTR), value :: p
    end subroutine free
  end interface
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

  subroutine phist_comm_create(comm_ptr, ierr) bind(C,name='phist_comm_create')
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


  subroutine phist_comm_delete(comm_ptr,ierr) bind(C,name='phist_comm_delete')
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

  subroutine phist_comm_get_rank(comm_ptr, rank, ierr) bind(C,name='phist_comm_get_rank')
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

  subroutine phist_comm_get_size(comm_ptr, nprocs, ierr) bind(C,name='phist_comm_get_size')
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
  ! helper functions to allocate aligned memory
  subroutine allocate_aligned_1(n, array, ierr)
    use, intrinsic :: iso_c_binding
    !------------------------------------------------------------
    integer(kind=C_SIZE_T),            intent(in)   :: n(1)
    real(kind=8), pointer, contiguous, intent(out)  :: array(:)
    integer,                           intent(out)  :: ierr
    !------------------------------------------------------------
    type(C_PTR) :: rawMem = C_NULL_PTR
    integer(kind=C_SIZE_T), parameter :: alignement = 64
    integer(kind=C_SIZE_T), parameter :: double = 8
    !------------------------------------------------------------

    ierr = posix_memalign(rawMem, alignement, product(n)*double)
    if ( ierr/=0 .or. .not. c_associated(rawMem) ) then
      write(*,*) 'Error allocating ', n*double/2**20, 'MB array!'
      ierr = -1
      return
    end if
    call c_f_pointer(rawMem, array, n)
  end subroutine allocate_aligned_1

  subroutine allocate_aligned_2(n, array, ierr)
    use, intrinsic :: iso_c_binding
    !------------------------------------------------------------
    integer(kind=C_SIZE_T),            intent(in)   :: n(2)
    real(kind=8), pointer, contiguous, intent(out)  :: array(:,:)
    integer,                           intent(out)  :: ierr
    !------------------------------------------------------------
    type(C_PTR) :: rawMem = C_NULL_PTR
    integer(kind=C_SIZE_T), parameter :: alignement = 64
    integer(kind=C_SIZE_T), parameter :: double = 8
    !------------------------------------------------------------

    ierr = posix_memalign(rawMem, alignement, product(n)*double)
    if ( ierr/=0 .or. .not. c_associated(rawMem) ) then
      write(*,*) 'Error allocating ', n*double/2**20, 'MB array!'
      ierr = -1
      return
    end if
    call c_f_pointer(rawMem, array, n)
  end subroutine allocate_aligned_2

  subroutine allocate_aligned_3(n, array, ierr)
    use, intrinsic :: iso_c_binding
    !------------------------------------------------------------
    integer(kind=C_SIZE_T),            intent(in)   :: n(3)
    real(kind=8), pointer, contiguous, intent(out)  :: array(:,:,:)
    integer,                           intent(out)  :: ierr
    !------------------------------------------------------------
    type(C_PTR) :: rawMem = C_NULL_PTR
    integer(kind=C_SIZE_T), parameter :: alignement = 64
    integer(kind=C_SIZE_T), parameter :: double = 8
    !------------------------------------------------------------

    ierr = posix_memalign(rawMem, alignement, product(n)*double)
    if ( ierr/=0 .or. .not. c_associated(rawMem) ) then
      write(*,*) 'Error allocating ', n*double/2**20, 'MB array!'
      ierr = -1
      return
    end if
    call c_f_pointer(rawMem, array, n)
  end subroutine allocate_aligned_3


  !================================================================================
  ! helper functions to deallocate aligned memory
  subroutine deallocate_aligned(array)
    use, intrinsic :: iso_c_binding, only: c_loc
    real(kind=8), target :: array(*)
    !--------------------------------------------------------------------------------

    call free(c_loc(array(1)))

  end subroutine deallocate_aligned


end module env_module

! just used to prevent compiler warnings about unused variables
subroutine TOUCH_EMPTY_FUNCTION(x)
end subroutine TOUCH_EMPTY_FUNCTION
