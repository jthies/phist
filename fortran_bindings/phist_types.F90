!! this module declares the basic KIND parameters
!! for arguments to PHIUST functions.
!!
!! the index types should be declared like this:
!!
!! INTEGER(lidx) :: nloc
!! INTEGER(gidx) :: nglob
!!
!! in order to also make use of the typical
!! PHIST data structure handles like TYPE(sdMat_ptr),
!! one should #include "phist_fort.h" and use the
!! macros defined there.
module phist_types

use, intrinsic :: iso_c_binding

include 'phist_enums.F90'

!! TODO: make this depend on the PHIST installation
integer, parameter :: lidx = c_int32_t
integer, parameter :: gidx = c_int64_t

end module phist_types
