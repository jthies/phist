extern "C" {

void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
  //TODO: maybe the halocommInit call could go here?
  *iflag=0;
  
  // TODO Christie compute row scaling and return it in *work as a vector or array
  
  // CARP sweep not implemented, we indicate this here already
  // to avoid confusion in the tests
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

//TODO Jonas, change interface to store real and imag part consecutively
void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
*iflag=PHIST_NOT_IMPLEMENTED;
return;
// TODO Christie, update this for X_i==NULL and ignoring sigmas
#if 0
    ghost_densemat_halo_comm_t comm = GHOST_DENSEMAT_HALO_COMM_INITIALIZER;
    PHIST_CHK_GERR(x->halocommInit(x,&comm),*iflag);
    PHIST_CHK_GERR(x->halocommStart(x,&comm),*iflag);
    PHIST_CHK_GERR(x->halocommFinalize(x,&comm),*iflag);
    
    PHIST_CHK_GERR(mat->kacz(mat,x,b,omega,1),*iflag);
    PHIST_CHK_GERR(x->averageHalo(x),*iflag);
    PHIST_CHK_GERR(mat->kacz(mat,x,b,omega,0),*iflag);
    PHIST_CHK_GERR(x->averageHalo(x),*iflag);

  *iflag=0;
  return;
#endif
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
  // can be ignored for now (related to JD)
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A,
void* work, int *iflag)
{
  // TODO Christie, delete the vector of row scaling elements
  *iflag=0;
  return;
}

} // extern "C"
