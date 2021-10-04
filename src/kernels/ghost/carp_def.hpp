/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include <ghost.h>
#include <ghost/helper.h>
#include <iostream>
#include <fstream>

extern "C" {

void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) vA, int numShifts, 
        _ST_ const sigma[],
        void** work, int* iflag)
{
// this only works with the GHOST kacz branch so far
  *iflag=PHIST_NOT_IMPLEMENTED;
#if 0
  //TODO: maybe the halocommInit call could go here?

#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  PHIST_CAST_PTR_FROM_VOID( ghost_sparsemat, A, vA, *iflag);

  
  double *rownorm = new double[A->nrows];
 
  for(int row=0; row < A->nrows; ++row) {
      rownorm[row]=0;
  }

  ghost_sell *sellmat = SELL(A);
  double *mval = (double *)sellmat->val;
  ghost_lidx idx;

  for(int row=0; row < A->nrows; ++row) {
      idx =  sellmat->chunkStart[row];
      for (int j=0; j<sellmat->rowLen[row]; ++j) {
            rownorm[row] += mval[idx]*mval[idx];
       }
   }

  //do dissection of zones
  ghost_rcm_dissect(A);

 *work = (void*)rownorm;

  
  // TODO Christie compute row scaling and return it in *work as a vector or array
    
  // CARP sweep not implemented, we indicate this here already
  // to avoid confusion in the tests
#endif
  return;
}

//TODO Jonas, change interface to store real and imag part consecutively
void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) vA,
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X,
        void* const work,
        _MT_ const * omega, int* iflag)
{

#ifdef PHIST_TESTING
  // check if the in/output and rhs have the correct maps. For Christies implementation of the kernel,
  // we need X and B to live in the *column map* of A, the actual ordering of the columns, not of an input
  // vector to the spMVM (the domain map).
  {
    phist_const_map_ptr map_x, map_b=NULL, col_map_A;
    PHIST_CHK_IERR(SUBR(sparseMat_get_col_map)(vA,&col_map_A,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_get_map)(X,&map_x,iflag),*iflag);
    // x and y must be correctly partitioned and permuted at this point, so demand *iflag=0 here:
    PHIST_CHK_IERR(phist_maps_compatible(map_x, col_map_A,iflag),*iflag);
    if (Rhs)
    {
      PHIST_CHK_IERR(SUBR(mvec_get_map)(Rhs,&map_b,iflag),*iflag);
      PHIST_CHK_IERR(phist_maps_compatible(map_b, col_map_A,iflag),*iflag);
    }
    
  }
#endif

// this only works with the GHOST kacz branch so far
  *iflag=PHIST_NOT_IMPLEMENTED;
#if 0
 *iflag=0;
 PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat, A, vA, *iflag);

 ghost_densemat *b = NULL;

if(Rhs != NULL)
 b=(ghost_densemat*)(Rhs);

 PHIST_CAST_PTR_FROM_VOID(ghost_densemat, x, X, *iflag);

 double omega_ = *omega;
 ghost_kacz_opts opts = GHOST_KACZ_OPTS_INITIALIZER;
 opts.omega = &omega_;
 opts.normalize = no;
 opts.direction = GHOST_KACZ_DIRECTION_FORWARD;
 ghost_kacz_rb(x,A,b,opts);
 opts.direction = GHOST_KACZ_DIRECTION_BACKWARD;
 ghost_kacz_rb(x,A,b,opts);
#endif
return;
// TODO Christie, update this for X_i==NULL and ignoring sigmas
/*#if 0
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
#endif*/
}

#ifndef IS_COMPLEX

// variant with real matrix and complex shifts is not implemented yet

void SUBR(carp_setup_rc)(TYPE(const_sparseMat_ptr) vA, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

void SUBR(carp_sweep_rc)(TYPE(const_sparseMat_ptr) vA,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#endif


void SUBR(carp_sweep_aug)(TYPE(const_sparseMat_ptr) A,
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) Q,
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X,
        TYPE(sdMat_ptr) q,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  // can be ignored for now (related to JD)
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

void SUBR(carp_sweep_aug_rc)(TYPE(const_sparseMat_ptr) A,
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
#if 0
  // TODO Christie, delete the vector of row scaling elements
   PHIST_ENTER_FCN(__FUNCTION__);
   *iflag=0;
   PHIST_CAST_PTR_FROM_VOID(double,rownorm,work,*iflag);
   delete[] rownorm;
#endif

  *iflag=0;
  return;
}

} // extern "C"
