// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D)). If we ever
// decide to add more data types (such as quad precision),
// this will help maintain existing implementations for the
// classical data types.

// get rid of all those warnings!
void SUBR(type_avail)(int *iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm,
        const char* filename,int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm,
        const char* filename,int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm,
        const char* filename,int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    const_map_ptr_t map, lidx_t nvec, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, const_comm_ptr_t comm, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  *iflag=-99;
}


void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_identity)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  *iflag=-99;
}


void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=-99;
}

void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=-99;
}


void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                          TYPE(const_sdMat_ptr) W, 
                              _ST_ beta,        TYPE(sdMat_ptr) C,
                              int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) M, int* rank, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  *iflag=-99;
}

#ifdef PHIST_KERNEL_LIB_BUILTIN
void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__)
  *iflag=-99;
}

void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__)
  *iflag=-99;
}
#endif

void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int* iflag)
{
  *iflag=-99;
}

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Dmvec_t* reV, Dmvec_t* imV, int *iflag)
{
  *iflag=-99;
}
# else
void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Smvec_t* reV, Smvec_t* imV, int *iflag)
{
  *iflag=-99;
}
# endif
#endif

