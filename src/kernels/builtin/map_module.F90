/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

!> \file map_module.f90
!! Defines map_module, the phist_map implementation of the builtin kernels
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!! \author "Jonas Thies <j.thies@tudelft.nl>
!!

#include "phist_config_fortran.h"
#include "phist_kernel_flags.h"
#include "phist_defs.h"

!> simple macro to check integer flags
#define CHECK_IFLAG(_flag,_value) (IAND(_flag,_value)==(_value))


!> Implements phist_map_* for builtin kernels
module map_module
  implicit none
  private

  public :: Map_t
  public :: map_setup
  public :: map_compatible_map
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
    !! offset array, length (0:nProcs)
    integer(kind=8), allocatable :: distrib(:)
    !! for each MPI process, contains the number of local elements (0:nProcs-1)
    integer,         allocatable :: nlocal(:)
    !! for global reorderings (1:nlocal(me)), not allocated by default
    integer(kind=8), allocatable :: global_idx(:)
    
    ! the may map contain coloring information:
    ! to access the elements of a vector (or rows of
    ! a matrix) one color at a time, do something like
    ! this:
    !
    ! do ic=1,map%nColors
    !   do jc=map%color_offset(ic),map%color_offset(ic+1)-1
    !     x(map%color_idx(jc))=...
    !
    ! This module just sets map%nColors to 0 and leaves the arrays unallocated
    ! but deletes them if they are allocated in map_delete().
    integer :: coloringType
    integer :: nColors
    integer,         allocatable :: color_offset(:), color_idx(:)
  end type Map_t

contains

  !================================================================================
  ! setup map from comm and number of rows
  subroutine map_setup(map, comm, n_glob, verbose, ierr)
    use mpi_f08
    !------------------------------------------------------------
    type(Map_t),    intent(inout) :: map
    integer,        intent(in)  :: comm
    integer(kind=8),intent(in)  :: n_glob
    logical,        intent(in)  :: verbose
    integer,        intent(out) :: ierr
    !------------------------------------------------------------
    integer :: i
    !------------------------------------------------------------

    map%comm = comm
    call mpi_comm_size(map%comm, map%nProcs, ierr)
    call mpi_Comm_rank(map%comm, map%me, ierr)

    if( verbose .and. map%me .eq. 0 ) then
      write(*,*) 'map%nProcs', map%nProcs
    end if

    allocate(map%distrib(0:map%nProcs))
    allocate(map%nlocal(0:map%nProcs-1))

    do i = 0, map%nProcs-1, 1
      map%nlocal(i) = int(n_glob/map%nProcs)
      if( i .lt. mod(n_glob,map%nProcs) ) then
        map%nlocal(i) = map%nlocal(i) + 1
      end if
    end do

    map%distrib(0) = 1
    do i = 1, map%nProcs, 1
      map%distrib(i) = map%distrib(i-1) + map%nlocal(i-1)
    end do

#ifdef PHIST_TESTING
if( verbose ) then
  write(*,*) map%me, 'setting up map with nglob = ', n_glob, ', nProcs = ', map%nProcs, &
    &  ', distrib = ', map%distrib, ', nlocal = ', map%nlocal
  flush(6)
end if
#endif
    if( map%distrib(map%nProcs) .ne. n_glob+1 ) then
      ierr = -1
    else
      ierr = 0
    end if

  map%nColors=0
  map%coloringType=0

  end subroutine map_setup


  !================================================================================
  ! check if two maps are compatible
  function map_compatible_map(map1, map2, reorder) result(res)
    use mpi
    !------------------------------------------------------------
    type(Map_t),  intent(in)            :: map1, map2
    logical,      intent(in), optional  :: reorder
    logical                             :: res
    !------------------------------------------------------------
    logical :: localRes
    integer :: ierr
    !------------------------------------------------------------

    res = .true.
    if( map1%comm   .ne. map2%comm   ) res = .false.
    if( map1%nProcs .ne. map2%nProcs ) res = .false.
    if( map1%me     .ne. map2%me     ) res = .false.

    if( .not. allocated(map1%nlocal) .or. .not. allocated(map2%nlocal) ) res = .false.

    if( .not. present(reorder) ) then

      ! just check local dimensions
      if( res ) then
        if( any(map1%nlocal .ne. map2%nlocal) ) res = .false.
      end if

    else

      ! check global number of elements
      if( map1%distrib(map1%nProcs) .ne. map2%distrib(map2%nProcs) ) res = .false.

      if( .not. reorder ) then
        ! no reordering allowed, must be identical!
        if( any(map1%distrib .ne. map2%distrib) ) res = .false.

        if( allocated(map1%global_idx) .and. allocated(map2%global_idx) ) then

          if( any(map1%global_idx .ne. map2%global_idx)) then
            res = .false.
          end if
          localRes = res
          call MPI_Allreduce(localRes, res, 1, MPI_LOGICAL, MPI_LAND, map1%comm, ierr);

        else if( allocated(map1%global_idx) .neqv. allocated(map2%global_idx) ) then
          res = .false.
        end if

      end if

    end if

  end function map_compatible_map


  !================================================================================
  ! some simple wrapper routines

  subroutine phist_map_create(map_ptr, comm_ptr, n_glob, ierr) bind(C,name='phist_map_create')
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),        intent(out) :: map_ptr
    type(C_PTR),        value       :: comm_ptr
    integer(C_INT64_T), value       :: n_glob
    integer(C_INT),     intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    integer, pointer :: comm
    logical :: verbose
    !------------------------------------------------------------
    verbose = .not. CHECK_IFLAG(ierr,PHIST_SPARSEMAT_QUIET)

    call c_f_pointer(comm_ptr,comm)

    allocate(map)

    call map_setup(map, comm, n_glob, verbose, ierr)

    map_ptr = c_loc(map)

  end subroutine phist_map_create


  subroutine phist_map_delete(map_ptr, ierr) bind(C,name='phist_map_delete')
    use, intrinsic :: iso_c_binding
    !------------------------------------------------------------
    type(C_PTR),    value       :: map_ptr
    integer(C_INT), intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    if (allocated(map%distrib)) then
      deallocate(map%distrib)
    end if
    if (allocated(map%nlocal)) then
      deallocate(map%nlocal)
    end if
    if (allocated(map%color_offset)) then
      deallocate(map%color_offset)
    end if
    if (allocated(map%color_idx)) then
      deallocate(map%color_idx)
    end if
    if (allocated(map%global_idx)) then
      deallocate(map%global_idx)
    end if
    deallocate(map)
    ierr = 0
  end subroutine phist_map_delete


  subroutine phist_map_get_comm(map_ptr, comm_ptr, ierr) bind(C,name='phist_map_get_comm')
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


  subroutine phist_map_get_local_length(map_ptr, nloc, ierr) bind(C,name='phist_map_get_local_length')
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

  subroutine phist_map_get_global_length(map_ptr, nglob, ierr) bind(C,name='phist_map_get_global_length')
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),        value       :: map_ptr
    integer(C_INT64_T), intent(out) :: nglob
    integer(C_INT),     intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    nglob = map%distrib(map%nProcs)-1
    ierr = 0
  end subroutine phist_map_get_global_length


  subroutine phist_map_get_ilower(map_ptr, ilower, ierr) bind(C,name='phist_map_get_ilower')
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),        value       :: map_ptr
    integer(C_INT64_T), intent(out) :: ilower
    integer(C_INT),     intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    ilower = map%distrib(map%me)-1
    ierr = 0
  end subroutine phist_map_get_ilower


  subroutine phist_map_get_iupper(map_ptr, iupper, ierr) bind(C,name='phist_map_get_iupper')
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),        value       :: map_ptr
    integer(C_INT64_T), intent(out) :: iupper
    integer(C_INT),     intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map
    !------------------------------------------------------------

    call c_f_pointer(map_ptr, map)
    iupper = map%distrib(map%me+1)-2
    ierr = 0
  end subroutine phist_map_get_iupper

  subroutine phist_maps_compatible(map1_ptr, map2_ptr, ierr) bind(C,name='phist_maps_compatible')
    use, intrinsic :: iso_c_binding
    use mpi
    !------------------------------------------------------------
    type(C_PTR),        value       :: map1_ptr, map2_ptr
    integer(C_INT),     intent(out) :: ierr
    !------------------------------------------------------------
    type(Map_t), pointer :: map1,map2
    !------------------------------------------------------------
    
    integer :: i, mpi_ierr, my_ierr, global_ierr
    logical :: map1_permuted, map2_permuted

    ! check if the two pointers refer to the same memory location
    if ( c_associated(map1_ptr,map2_ptr) ) then
      ierr=0
      return
    end if

    call c_f_pointer(map1_ptr, map1)
    call c_f_pointer(map2_ptr, map2)
    
    if (.not. (associated(map1) .and. associated(map2))) then
      ierr=PHIST_BAD_CAST
      return
    end if
    
    ierr=-1
    
    ! check if MPI setup is compatible
    if ( map1%comm      /= map2%comm    .or. &
         map1%nprocs   /= map2%nprocs  .or. &
         map1%me       /= map2%me ) then
      return
    end if
    
    ! check global size is the same
    if ( map1%distrib(map1%nprocs) /= map2%distrib(map2%nprocs) ) then
      return
    end if

    map1_permuted=allocated(map1%global_idx)
    map2_permuted=allocated(map2%global_idx)

    ! check compatibility if both or none of the maps have been reordered
    if ( map1_permuted .eqv. map2_permuted ) then
        
      do i=0,map1%nProcs
        if (map1%distrib(i)/=map2%distrib(i) ) then
          return
        end if
      end do
      ! maps are the same if both maps are unpermuted
      if (.not. map1_permuted) then
        ierr=0
        return
      end if
    end if
    
    if (map1_permuted .and. map2_permuted) then
      ! check if the permutation is the same
      ierr=0
      do i=1,map1%nlocal(map1%me)
        if ( map1%global_idx(i)/=map2%global_idx(i) ) then
          ierr=1
          exit
        end if
      end do
      my_ierr=ierr
      call MPI_Allreduce(my_ierr,global_ierr,1,MPI_INTEGER,MPI_MAX,map1%comm,mpi_ierr)
      if (mpi_ierr/=MPI_SUCCESS) then
        ierr=PHIST_MPI_ERROR
        return
      end if
      ierr=global_ierr
      if (ierr/=0) then
        ! we don't support this type of construction (both maps permuted with different orderings)
        ierr=-1
      end if
      return
    end if
    
    ! this is the standard case, one map is the permutation of the other.
    ! we do not distinguish between global and local, though. Note that we
    ! already verified that the maps have the same global size above.
    ierr=2
    
  end subroutine phist_maps_compatible

  


end module map_module
