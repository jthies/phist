!! this header defines the basic data types that appear as 
!! arguments to PHIUST functions.
!!
!! the interface objects should be declared this
!!
!! TYPE(comm_ptr)  :: comm
!! TYPE(DsdMat_ptr) :: M
!! TYPE(Zconst_mvec_ptr) :: V
!! etc.
!!
!! in fact, the latter all evaluate to c_ptr, so it does not matter too much
!! which declaration is used, but the code is more readable if the declaration
!! contains information about const-ness and scalar data type of the target object.
!!
!! For clarity of the generated interfaces, the type specifiers are added, as in 
!! the example TYPE(DsdMat_ptr) above. We also provide a macro for each type without
!! the prefix, which allows the more phist-like syntax TYPE(sdMat_ptr) etc. when 
!!declaring objects
!!
#ifndef PHIST_FORT_H
#define PHIST_FORT_H

#include "phist_config.h"
#include "phist_kernel_flags.h"

!! eqivalen of phist_void_aliases.h: since there are no
!! typedefs in Fortran, we need to use macros
#define comm_ptr c_ptr
#define const_comm_ptr c_ptr
#define map_ptr c_ptr
#define const_map_ptr c_ptr
#define context_ptr c_ptr
#define const_context_ptr c_ptr

#define sdMat_ptr  c_ptr
#define SsdMat_ptr  c_ptr
#define CsdMat_ptr  c_ptr
#define DsdMat_ptr  c_ptr
#define ZsdMat_ptr  c_ptr

#define const_sdMat_ptr  c_ptr
#define Sconst_sdMat_ptr  c_ptr
#define Cconst_sdMat_ptr  c_ptr
#define Dconst_sdMat_ptr  c_ptr
#define Zconst_sdMat_ptr  c_ptr

#define mvec_ptr  c_ptr
#define Smvec_ptr  c_ptr
#define Cmvec_ptr  c_ptr
#define Dmvec_ptr  c_ptr
#define Zmvec_ptr  c_ptr

#define const_mvec_ptr  c_ptr
#define Sconst_mvec_ptr  c_ptr
#define Cconst_mvec_ptr  c_ptr
#define Dconst_mvec_ptr  c_ptr
#define Zconst_mvec_ptr  c_ptr

#define sparseMat_ptr  c_ptr
#define SsparseMat_ptr  c_ptr
#define CsparseMat_ptr  c_ptr
#define DsparseMat_ptr  c_ptr
#define ZsparseMat_ptr  c_ptr

#define const_sparseMat_ptr  c_ptr
#define Sconst_sparseMat_ptr  c_ptr
#define Cconst_sparseMat_ptr  c_ptr
#define Dconst_sparseMat_ptr  c_ptr
#define Zconst_sparseMat_ptr  c_ptr

#endif
