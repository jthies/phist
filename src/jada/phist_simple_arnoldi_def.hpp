//! a simple Arnoldi process to start up the JaDa iteration.
//! Given a minimum basis size m, compute V(:,1:m), 
//! H(1:m+1,1:m) such that A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
//! input: v0 with ||v0||_2=1, V and H allocated with m+1 columns
//! and nloc resp. m+1 rows.
//!
//! TODO - block Arnoldi, cf. matlab/krylov/arnoldi.m for a prototype.
//! TODO - we may want to include a check if any Ritz values have already
//!        converged in Arnoldi, but this requires some additional programming
//!        effort which we leave for later.
void SUBR(simple_arnoldi)(TYPE(const_op_ptr) op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(sdMat_ptr) H, int m, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  TYPE(mvec_ptr) v=NULL,av=NULL,vprev=NULL;
  TYPE(sdMat_ptr) R1=NULL,R2=NULL;
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V,v0,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&v,0,0,ierr),*ierr);
  _MT_ v0norm;
  PHIST_CHK_IERR(SUBR(mvec_normalize)(v,&v0norm,ierr),*ierr);
  
  for (int i=0;i<m-1;i++)
    {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&vprev,0,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&av,i+1,i+1,ierr),*ierr);
    PHIST_CHK_IERR(op->apply(st::one(),op->A,v,st::zero(),av,ierr),*ierr);
    // orthogonalize, Q*R1 = W - V*R2
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R2,0,i,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,i+1,i+1,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(orthog)(vprev,av,R1,R2,2,ierr),*ierr);
    v=av; av=NULL;
    }
  }


