void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}


void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

void SUBR(carp_sweep_aug)(TYPE(const_sparseMat_ptr) A,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Q,
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        TYPE(sdMat_ptr) q_r, TYPE(sdMat_ptr) q_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A,
        void* work, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}



