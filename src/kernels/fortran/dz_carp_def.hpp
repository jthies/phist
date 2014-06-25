#ifndef IS_DOUBLE
#error "this file is only intended for the Z-case"
#endif
#ifndef IS_COMPLEX
#error "this file is only intended for the Z-case"
#endif

void SUBR(carp_setup)(TYPE(const_crsMat_ptr) A, _ST_ const sigma[], 
        Dmvec_t *nrms, TYPE(mvec_ptr)* importVec, int* ierr)
{
  *ierr=-99;
  return;
}

void SUBR(dkswp)(TYPE(const_crsMat_ptr) A, _ST_ const sigma[], 
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X, 
        Dmvec_t const*  nrm_ai2, TYPE(mvec_ptr) importVec,
        _MT_ const * omega, int* ierr)
{
  *ierr=-99;
  return;
}

// mixed variants: real matrix with complex shifts

void phist_DZcarp_setup(DcrsMat_t const*  A, _ST_ const sigma[], 
        Dmvec_t* nrms, Zmvec_t** importVec, int* ierr)
{
  //TODO
  *ierr=-99;
  return;
}

void phist_DZdkswp(DcrsMat_t const* A, _ST_ const sigma[], 
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X, 
        Dmvec_t const* nrm_ai2, Zmvec_t* importVec,
        _MT_ const * omega, int* ierr)
{
  //TODO
  *ierr=-99;
  return;
}

