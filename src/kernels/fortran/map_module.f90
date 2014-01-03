module map_module
  implicit none
  private

  public :: Map_t
  public :: map_setup
  !public :: phist_map_create
  !public :: phist_map_delete
  !public :: phist_map_get_comm
  !public :: phist_map_get_local_length
  !public :: phist_map_get_ilower
  !public :: phist_map_get_iupper


  type Map_t
    integer :: comm
    integer :: nProcs
    integer :: me
    integer(kind=8), allocatable :: distrib(:)
    integer,         allocatable :: nlocal(:)
  end type Map_t

contains

  !================================================================================
  ! setup map from comm and number of rows
  subroutine map_setup(map, comm, n_glob, ierr)
    use mpi
    !------------------------------------------------------------
    type(Map_t),    intent(out) :: map
    integer,        intent(in)  :: comm
    integer(kind=8),intent(in)  :: n_glob
    integer,        intent(out) :: ierr
    !------------------------------------------------------------
    integer :: i
    !------------------------------------------------------------

    map%comm = comm
    call mpi_comm_size(map%comm, map%nProcs, ierr)
    call mpi_Comm_rank(map%comm, map%me, ierr)

    allocate(map%distrib(0:map%nProcs))
    allocate(map%nlocal(0:map%nProcs-1))

    do i = 0, map%nProcs-1
      map%nlocal(i) = int(n_glob/map%nProcs)
      if( i .lt. mod(n_glob,map%nProcs) ) then
        map%nlocal(i) = map%nlocal(i) + 1
      end if
    end do

    map%distrib(0) = 1
    do i = 1, map%nProcs
      map%distrib(i) = map%distrib(i-1) + map%nlocal(i-1)
    end do

    if( map%distrib(map%nProcs) .ne. n_glob+1 ) then
      ierr = -1
    else
      ierr = 0
    end if

  end subroutine map_setup


  !================================================================================
  ! some simple wrapper routines

  subroutine phist_map_create(map_ptr, comm_ptr, n_glob, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),        intent(out) :: map_ptr
    type(C_PTR),        value       :: comm_ptr
    integer(C_INT64_T), value       :: n_glob
    integer(C_INT),     intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    allocate(map)

    call map_setup(map, MPI_COMM_WORLD, n_glob, ierr)

    map_ptr = c_loc(map)

  end subroutine phist_map_create


  subroutine phist_map_delete(map_ptr, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    !------------------------------------------------------------
    type(C_PTR),    value       :: map_ptr
    integer(C_INT), intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    deallocate(map)
    ierr = 0
  end subroutine phist_map_delete


  subroutine phist_map_get_comm(map_ptr, comm_ptr, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    value       :: map_ptr
    type(C_PTR),    intent(out) :: comm_ptr
    integer(C_INT), intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    comm_ptr = c_loc(map%comm)
    ierr = 0
  end subroutine phist_map_get_comm


  subroutine phist_map_get_local_length(map_ptr, nloc, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    value       :: map_ptr
    integer(C_INT), intent(out) :: nloc, ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    nloc = map%nlocal(map%me)
    ierr = 0
  end subroutine phist_map_get_local_length


  subroutine phist_map_get_ilower(map_ptr, ilower, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    value       :: map_ptr
    integer(C_INT), intent(out) :: ilower, ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    ilower = map%distrib(map%me)-1
    ierr = 0
  end subroutine phist_map_get_ilower


  subroutine phist_map_get_iupper(map_ptr, iupper, ierr) bind(C)
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),    value       :: map_ptr
    integer(C_INT), intent(out) :: iupper, ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    iupper = int(map%distrib(map%me+1)-1)
    ierr = 0
  end subroutine phist_map_get_iupper


end module map_module
