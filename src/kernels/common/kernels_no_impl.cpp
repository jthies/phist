// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D)). If we ever
// decide to add more data types (such as quad precision),
// this will help maintain existing implementations for the
// classical data types.
//
// This file disables all kernel functions for the given data type,
// it is also possible to just include selected kernels_no_* files 
// if some features are not implemented. The other files typically 
// provide default implementations that will work but not be optimal
// in any way.

// get rid of all those warnings!
extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, phist_const_comm_ptr comm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
        phist_sparseMat_rowFunc rowFunPtr, void* last_arg, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    phist_const_map_ptr map, phist_lidx nvec, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, phist_const_map_ptr map, 
    _ST_* values, phist_lidx lda, int nvec,
    int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, phist_const_comm_ptr comm, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, phist_const_comm_ptr comm,
        _ST_* values, phist_lidx lda, int nrows, int ncols,
        int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, phist_const_map_ptr* map, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#ifdef PHIST_HIGH_PRECISION_KERNELS
extern "C" void SUBR(sdMat_extract_error)(TYPE(sdMat_ptr) M, _ST_** MC_raw, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
        phist_mvec_elemFunc funPtr, void* last_arg, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif
extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr) x, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                          TYPE(const_sdMat_ptr) W, 
                              _ST_ beta,        TYPE(sdMat_ptr) C,
                              int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#ifdef PHIST_KERNEL_LIB_BUILTIN
extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__)
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__)
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Dmvec* reV, phist_Dmvec* imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Smvec* reV, phist_Smvec* imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
# endif
#endif

#include "kernels_no_inplace_VC.cpp"
#include "kernels_no_VC_add_WD.cpp"
#include "kernels_no_carp.cpp"
#include "kernels_no_gpu.cpp"
#include "kernels_no_fused.cpp"
