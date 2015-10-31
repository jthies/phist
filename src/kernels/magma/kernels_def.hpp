/*! \file magma/kernels_def.hpp
 * included by magma/kernels.cpp
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/


extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    const_map_ptr_t map, lidx_t nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, const_comm_ptr_t comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}
                  

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#ifdef PHIST_HIGH_PRECISION_KERNELS
extern "C" void SUBR(sdMat_extract_error)(TYPE(sdMat_ptr) V, _ST_** err, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif
extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
        int (*funPtr)(ghost_gidx_t,ghost_lidx_t,void*), int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_times_sdMat_augmented)(_ST_ alpha,  TYPE(const_mvec_ptr)  V, 
                                                              TYPE(const_sdMat_ptr) C, 
                                                 _ST_ beta,   TYPE(mvec_ptr)        W,
                                                              TYPE(sdMat_ptr)       D,
                                                 int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V, 
                                                            TYPE(const_sdMat_ptr) C,
                                                            TYPE(mvec_ptr) W, 
                                                            TYPE(const_sdMat_ptr) D,
                                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(mvecT_times_mvec_times_sdMat_inplace)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                                       TYPE(mvec_ptr)        W,
                                                                       TYPE(const_sdMat_ptr) C,
                                                           _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                           int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
}

// NOTE: see the description of sparseMat_read_mm on how we treat input flags for this function
extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne, void* last_arg,
                int (*rowFunPtr)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*,void*), 
                int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#include "../common/kernels_nogpu.c"


