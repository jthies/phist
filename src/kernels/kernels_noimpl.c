// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D)). If we ever
// decide to add more data types (such as quad precision),
// this will help maintain existing implementations for the
// classical data types.

// get rid of all those warnings!
void SUBR(type_avail)(int *ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* A, const_comm_ptr_t comm,
        const char* filename,int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* A, const_comm_ptr_t comm,
        const char* filename,int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* A, const_comm_ptr_t comm,
        const char* filename,int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    const_map_ptr_t map, lidx_t nvec, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, const_comm_ptr_t comm, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* ierr)
{
  *ierr=-99;
}


void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_view_scattered)(TYPE(mvec_ptr) V, TYPE(mvec_ptr)* Vscat,
    int* cols, int ncols, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_delete)(TYPE(crsMat_ptr) A, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* ierr)
{
  *ierr=-99;
}


void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMatT_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr)
{
  *ierr=-99;
}

void SUBR(crsMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr)
{
  *ierr=-99;
}


void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
{
  *ierr=-99;
}

void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
{
  *ierr=-99;
}

void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* ierr)
{
  *ierr=-99;
}

#ifdef PHIST_KERNEL_LIB_FORTRAN
void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *ierr)
{
  *ierr=-99;
}

void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *ierr)
{
  *ierr=-99;
}
#endif

void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int* ierr)
{
  *ierr=-99;
}

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Dmvec_t* reV, Dmvec_t* imV, int *ierr)
{
  *ierr=-99;
}
# else
void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Smvec_t* reV, Smvec_t* imV, int *ierr)
{
  *ierr=-99;
}
# endif
#endif

