//                                                               
// This is the main orthogonalization routine in PHIST.          
// It takes an orthogonal basis V and a set of vectors W,        
// and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.        
// (Q overwrites W here).                                        
// The matrices R1 and R2 must be pre-allocated by the caller.   
//                                                               
// The algorithm used is up to numSweeps steps of classical      
// block Gram-Schmidt, alternated with QR factorizations to      
// normalize W. The method stops if the reduction in norm of W   
// by a CGS step is less than approx. a factor sqrt(2).          
//                                                               
// If we find that W does not have full column rank,             
// the matrix Q is augmented with random vectors which are made  
// mutually orthogonal and orthogonal against V. The original    
// rank of W is returned in ierr at the end of the routine. If   
// it happens somewhere during the process, we return an error   
// code ierr=-7.                                                 
//                                                               
// If a breakdown occurs, indicating that one of the columns of  
// W lives in the space spanned by the columns of V, ierr=-8 is  
// returned. A more convenient behavior may be added later, like 
// randomizing the column(s) as before.                          
//                                                               
// If the decrease in norm in one of the columns                 
// in the last CGS sweep indicates that the algorithm has not    
// yet converged, we return ierr=-9 to indicate that more steps  
// may be advisable. This should not happen in practice if       
// numSweeps>=2 ('twice is enough').                             
//                                                               
void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps,
                     int* ierr)
  {
#include "phist_std_typedefs.hpp"

  int m,k,i,j;
  MT* normW0, *normW1;
  MT breakdown; // to check for V'W becoming near 0 for some vector in W
  int rankW;
  bool stopGS=false;
  
  // auxiliary matrices
  st::sdMat_t *R1p,*R2p,*R1pp;
  const_comm_ptr_t comm;

  *ierr=0;
  
  // all vectors and matrices must be allocated.
  if (V==NULL || W==NULL || R1==NULL || R2==NULL) 
    {
    *ierr=-1;
    return;
    }

  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,ierr),*ierr);

  if (k==0) // no vectors to be orthogonalized
    {
    return;
    }

  if (m==0)
    {
    PHIST_CHK_IERR(SUBR(mvec_QR)(W,R1,ierr),*ierr);
    return;
    }

  if (numSweeps>1)
    {
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1p,k,k,comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1pp,k,k,comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R2p,m,k,comm,ierr),*ierr);
    }

#ifdef TESTING
  // check that all array dimensions are correct
  PHIST_CHK_IERR((int)(V==NULL || W==NULL || R1==NULL || R2==NULL),*ierr);
  lidx_t n,tmp;
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(W,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((n==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R1,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((k==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R1,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((k==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R2,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((m==tmp)?0:-1),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R2,&tmp,ierr),*ierr);
  PHIST_CHK_IERR(((k==tmp)?0:-1),*ierr);
#endif

  // compute the norms of the columns in W (for checking the result later)
  normW0=new MT[k];
  normW1=new MT[k];
  
  // determine original norms of W vectors and breakdown tolerance
  PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW0,ierr),*ierr);
  breakdown = normW0[0];
  for (i=1;i<k;i++) std::min(normW0[i], breakdown);
  breakdown*=mt::eps();
  
  // orthogonalize against V (first CGS sweep)

  //R2=V'*W;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2,ierr),*ierr);
  //W=W-V*R2;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,R2,st::one(),W,ierr),*ierr);
  // norms after first CGS sweep. This could be done cheaper by using R1 computed below,
  // since ||W_j||_2 = ||R1_j||_2 (TODO - that's how it's done in Belos).
  PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW1,ierr),*ierr);

  // orthogonalize W after the first GS pass. This gives us the rank  
  // of W (ierr>0 indicates the dimension of the null space of W, i.e.
  // the rank of W is k-ierr). If W is not full rank, it is augmented 
  // by random vectors. These still have to be orthogonalized against 
  // V, which we do next. The projection coefficients for this part   
  // are thrown away as the randomized vectors are not really related 
  //to W anyway.
  SUBR(mvec_QR)(W,R1,ierr);
  if (*ierr<0)
    {
    PHIST_OUT(0,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
        *ierr,(phist_retcode2str(*ierr)),"Xmvec_QR(W,R1,ierr)",__FILE__,__LINE__); 
    return;
    }
  if (*ierr>0) rankW=k-*ierr;
  if (rankW<k)
    {
    st::mvec_t *Wrnd=NULL;
    st::sdMat_t *Rrnd;
    int n0=*ierr;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(W,&Wrnd,rankW,k-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&Rrnd,m,n0,comm,ierr),*ierr);
    //R2=V'*W;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,Wrnd,st::zero(),Rrnd,ierr),*ierr);
    //W=W-V*R2;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,Rrnd,st::one(),Wrnd,ierr),*ierr);
    // throw away projection coefficients.
    PHIST_CHK_IERR(SUBR(sdMat_delete)(Rrnd,ierr),*ierr);
    }

  for (i=1;i<numSweeps;i++)
    {
    MT maxRed=1.0;
    for (j=0;j<k;j++)
      {
      maxRed=std::min(maxRed,normW1[j]/normW0[j]);
      normW0[j]=normW1[j];
      }

    if (maxRed>0.7)
      {
      stopGS=true;
      break;
      }

    //R2p=V'*W;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2p,ierr),*ierr);

    //W=W-V*R2';
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,R2p,st::one(),W,ierr),*ierr);

    //R2=R2+R2';
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R2p,st::one(),R2,ierr),*ierr);

    // orthogonalize W. The situation where the rank of W becomes
    // smaller than k here is (I think) very unlikely because we 
    // added random vectors before if rank(W) was not full. If   
    // it happens, we return with error code -10.                
    SUBR(mvec_QR)(W,R1p,ierr);
    if (*ierr<0)
      {
      PHIST_OUT(0,"Error code %d (%s) returned from call %s\n(file %s, line %d)",\
          *ierr,(phist_retcode2str(*ierr)),"Xmvec_QR(W,R1,ierr)",__FILE__,__LINE__); 
      return;
      }
    else if (*ierr>0)
      {
      PHIST_OUT(0,"Unexpected rank deficiency in orthog routine\n(file %s, line %d)",\
                __FILE__,__LINE__); 
      *ierr=-10;
      return;
      }
    
    //R1pp=R1
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1pp,ierr),*ierr);
    //R1=R1p*R1pp;
    PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R1p,R1pp,st::zero(),R1,ierr),*ierr);
    // again, the norm can be computed from the coefficients in R2p (TODO)
    PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW1,ierr),*ierr);    
    }

  // if in the last CGS sweep a column decreased too much in norm, 
  // return an error (we could return a warning, but since we used 
  // all positive numbers for the case of rank deficient W already,
  // the caller will have to deal with the -5 case specially.      
  if (stopGS==false)
    {
    // we need to check again, for instance, the 2nd pass may have given
    // the desired small reduction but the check would have been performed
    // at the beginning of pass 3, which was not done.
    for (i=0;i<k;i++)
      {
      if (normW1[i]<0.7*normW0[i])
        {
        // more orthogonalization steps may be required.
        *ierr=-9;
        }
      }
    }

  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1p,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1pp,ierr),*ierr);
  delete [] normW0;
  delete [] normW1;
  }
