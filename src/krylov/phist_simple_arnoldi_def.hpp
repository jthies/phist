//! a simple Arnoldi process to start up the JaDa iteration.
//! Given a minimum basis size m, compute V(:,1:m+1), 
//! H(1:m+1,1:m) such that A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
//! input: v0, V and H allocated with m+1 resp. m columns
//! and nloc resp. m+1 rows.
//!
//! TODO - block Arnoldi, cf. matlab/krylov/arnoldi.m for a prototype.
//! TODO - we may want to include a check if any Ritz values have already
//!        converged in Arnoldi, but this requires some additional programming
//!        effort which we leave for later.
void SUBR(simple_arnoldi)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, int m, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!");
    PHIST_CHK_IERR(*ierr = 99, *ierr);
  }
  TYPE(mvec_ptr) v=NULL,av=NULL,vprev=NULL;
  TYPE(sdMat_ptr) R1=NULL,R2=NULL;
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V,v0,0,0,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&v,0,0,ierr),*ierr);
  _MT_ v0norm;
  PHIST_CHK_IERR(SUBR(mvec_normalize)(v,&v0norm,ierr),*ierr);
  
  int nrH,ncH,ncV;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&ncV,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(H,&nrH,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(H,&ncH,ierr),*ierr);
  // note: we actually don't care about the dimensions of V and H as long as they are large 
  // enough since we just work with views in this function. However, we issue a warning 
  // (positive ierr) if the array dimensions are larger than expected so that the user 
  // doesn't unexpectedly end up with empty rows/cols
  if (ncV<m+1 || nrH<m+1 || ncH<m) {*ierr = -1; return;}
  if (ncV!=m+1 || nrH!=m+1 || ncH!=m)
    {
    PHIST_SOUT(PHIST_VERBOSE,"REMARK: input vectors/matrix to arnoldi are larger than necessary.");
    PHIST_SOUT(PHIST_VERBOSE,"        Requested subspace dimension is %d,",m);
    PHIST_SOUT(PHIST_VERBOSE,"        H is %dx%d (expecting %dx%d), V with %d cols (expecting %d)",
                                      nrH, ncH, m+1,m,ncV,m+1);
    }
  
  for (int i=0;i<m;i++)
    {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&vprev,0,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&av,i+1,i+1,ierr),*ierr);
    PHIST_CHK_IERR(A_op->apply(st::one(),A_op->A,v,st::zero(),av,ierr),*ierr);
    // orthogonalize, Q*R1 = W - V*R2
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R2,0,i,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,i+1,i+1,i,i,ierr),*ierr);
    PHIST_CHK_NEG_IERR(SUBR(orthog)(vprev,av,R1,R2,3,ierr),*ierr);
    if( *ierr > 0 )
    {
      PHIST_SOUT(PHIST_INFO,"found invariant subspace in arnoldi, expanding basis with a randomly generated orthogonal vector");
    }
    std::swap(v,av);     // swap the vectors v and av,av will be deleted and re-build as a 
                         // new view in the next step
    }
  // delete the views that we have created
  PHIST_CHK_IERR(SUBR(mvec_delete)(v,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(av,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(vprev,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R2,ierr),*ierr);
  }


