#include "phist_carp_decl.h"

extern "C" {

void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
  //TODO: maybe the halocommInit call could go here?
  *iflag=0;
  return;
}


void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
*iflag=-99;
return;
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

void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A,
void* work, int *iflag)
{
  *iflag=0;
  return;
}






} // extern "C"
