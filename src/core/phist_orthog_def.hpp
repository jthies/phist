// This is the main orthogonalization routine in PHIST.          
// It takes an orthogonal basis V and a set of vectors W,        
// and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.        
// (Q overwrites W here).                                        
// The matrices R1 and R2 must be pre-allocated by the caller.   
// If relax=0, we guarantee orthogonality to machine precision.  
// If not, we just check and report an orthogonalization failure 
// as ierr=1.                                                    
void _SUBR_(orthog)(_TYPE_(const_mvec_ptr) V,
                     _TYPE_(mvec_ptr) W,
                     _TYPE_(sdMat_ptr) R1,
                     _TYPE_(sdMat_ptr) R2,
                     int relax, int* ierr)
  {
#include "phist_std_typedefs.hpp"

  int m,k;
  PHIST_CHK_IERR(_SUBR_(mvec_num_vectors)(V,&m,ierr),*ierr);
  PHIST_CHK_IERR(_SUBR_(mvec_num_vectors)(W,&k,ierr),*ierr);

#ifdef TESTING
  // check that all pointers are not NULL and all array dimensions are correct
  PHIST_CHK_IERR((int)(V==NULL || W==NULL || R1==NULL || R2==NULL),*ierr);
  int n,tmp;
  PHIST_CHK_IERR(_SUBR_(mvec_my_length)(V,&n,ierr),*ierr);
  PHIST_CHK_IERR(_SUBR_(mvec_my_length)(W,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((n==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(_SUBR_(sdMat_get_nrows)(R1,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((k==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(_SUBR_(sdMat_get_ncols)(R1,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((k==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(_SUBR_(sdMat_get_nrows)(R2,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((m==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(_SUBR_(sdMat_get_ncols)(R2,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((k==tmp)?0:-1),*ierr);
#endif

  // orthogonalize against V

  //R2=V'*W;
  PHIST_CHK_IERR(_SUBR_(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2,ierr),*ierr);
  //W=W-V*R2;
  PHIST_CHK_IERR(_SUBR_(mvec_times_sdMat)(-st::one(),V,R2,st::one(),W,ierr),*ierr);

  if (relax==0)
    {
    // create temporary matrix R2'
    sdMat_ptr_t R2p;
    const_comm_ptr_t comm;
    PHIST_CHK_IERR(_SUBR_(mvec_get_comm)(V,&comm,ierr),*ierr);
    PHIST_CHK_IERR(_SUBR_(sdMat_create)(&R2p,m,k,comm,ierr),*ierr);
    //R2p=V'*W;
    PHIST_CHK_IERR(_SUBR_(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2p,ierr),*ierr);
    //W=W-V*R2';
    PHIST_CHK_IERR(_SUBR_(mvec_times_sdMat)(-st::one(),V,R2p,st::one(),W,ierr),*ierr);
    //R2=R2+R2';
    PHIST_CHK_IERR(_SUBR_(sdMat_add_sdMat)(st::one(),R2p,st::one(),R2,ierr),*ierr);

    PHIST_CHK_IERR(_SUBR_(sdMat_delete)(R2p,ierr),*ierr);
    }

  // orthogonalize W
  PHIST_CHK_IERR(_SUBR_(mvec_QR)(W,R1,ierr),*ierr);
  }
