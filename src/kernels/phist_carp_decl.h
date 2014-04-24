#ifdef __cplusplus
extern "C" {
#endif

//! create data structures needed for subsequent calls to dkswp: 
//! (1) compute the 2-norm (squared) of each row of sparse matrix A-sigma[j]*I. The
//! number of shifts sigma is determined by the number of columns in nrms, so nrms 
//! must be allocated beforehand by the user.
//! (2) importVec - an import vector. This is allocated internally.
//! Both nrms and importVec must be deleted by the user when no more calls to dkswp
//! are desired. If *importVec is not NULL on input, it is deleted and reallocated.
void SUBR(carp_setup)(TYPE(const_crsMat_ptr) A, _ST_ const sigma[], 
        Dmvec_t *nrms, TYPE(mvec_ptr)* importVec, int* ierr);

//! forward/backward sweep of Kaczmarz/CARP algorithm (SSOR sweep on the normal equations),
//! with matrix A-sigma[j]*I applied to vector column j. The number of systems is determined 
//! by the number of columns in  X, and all other in/out args must have compatible number of 
//! cols/entries.
void SUBR(dkswp)(TYPE(const_crsMat_ptr) A, _ST_ const sigma[], 
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X, 
        Dmvec_t const*  nrm_ai2, TYPE(mvec_ptr) importVec,
        _MT_ const * omega, int* ierr);

// mixed variants: real matrix with complex shifts
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
void phist_DZcarp_setup(DcrsMat_t const*  A, _ST_ const sigma[], 
        Dmvec_t* nrms, Zmvec_t** importVec, int* ierr);

void phist_DZdkswp(DcrsMat_t const* A, _ST_ const sigma[], 
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X, 
        Dmvec_t const* nrm_ai2, Zmvec_t* importVec,
        _MT_ const * omega, int* ierr);
# else
void phist_SCcarp_setup(ScrsMat_t const*  A, _ST_ const sigma[], 
        Smvec_t* nrms, Cmvec_t** importVec, int* ierr);

void phist_SCdkswp(ScrsMat_t const* A, _ST_ const sigma[], 
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X, 
        Smvec_t const* nrm_ai2, Cmvec_t* importVec,
        _MT_ const * omega, int* ierr);
# endif
#endif

#ifdef __cplusplus
} //extern "C"
#endif


