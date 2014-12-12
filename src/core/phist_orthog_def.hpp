// helper routine in case the kernel lib does not provide mvec_QR, which
// is quite an advanced kernel to ask for.
void SUBR(mvec_QR_fallback)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, _MT_* nrmsV, int *ierr);

//! orthogonalize an mvec against an already orthogonal one. \ingroup core

//! This is the main orthogonalization routine in PHIST.           
//! It takes an orthogonal basis V and a set of vectors W,         
//! and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.         
//! (Q overwrites W here).                                         
//! The matrices R1 and R2 must be pre-allocated by the caller.    
//!                                                                
//! The algorithm used is up to numSweeps steps of classical       
//! block Gram-Schmidt, alternated with QR factorizations to       
//! normalize W. The method stops if the reduction in norm of W    
//! by a CGS step is less than approx. a factor sqrt(2).           
//!                                                                
//! If we find that W-V*R2 does not have full column rank,         
//! the matrix Q is augmented with random vectors which are made   
//! mutually orthogonal and orthogonal against V. In this case the 
//! dimension of the null space of W-V*R2 is returned in ierr>0.   
//!                                                                
//! The implementation makes use of the kernel function mvec_QR    
//! for the inner orthogonalization of W. If this function is not  
//! implemented by the kernel lib (i.e. returns -99), we use a     
//! fallback variant based on the SVQB algorithm. In this case,    
//! the relation Q*R1 = W-V*R2 does not hold, but Q is orthonormal 
//! and orthogonal against V anyway with the nullspace replaced    
//! by random vectors. To indicate the invalid R's, ierr=-7 is     
//! returned (which means that the rank of V is lost to the out-   
//! side world, even though it was computed and the corresponding  
//! columns of Q were randomized).                                 
//!                                                                
//! If no random orthogonal vectors can be generated (after some   
//! tries) ierr=-8 is returned. This may indicate a problem with   
//! the random vector generator.                                   
//!                                                                
//! If the decrease in norm in one of the columns                  
//! in the last CGS sweep indicates that the algorithm has not     
//! yet converged, we return ierr=-9 to indicate that more steps   
//! may be advisable. This should not happen in practice if        
//! numSweeps>=2 ('twice is enough').                              
//!                                                                
void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps,
                     int* ierr)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  int m,k;
  MT breakdown; // to check for V'W becoming near 0 for some vector in W
  int rankW;
  bool stopGS=false;
  
  // auxiliary matrices
  st::sdMat_t *R1p,*R2p,*R1pp;

  const_comm_ptr_t comm=NULL;

  bool useSVQB=false;
  *ierr=0;
  
  if (V==NULL)
  {
    m=0;
  }
  else
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,ierr),*ierr);
  }
  if (W==NULL)
  {
    k=0;
  }
  else
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,ierr),*ierr);
  }
  if (k==0) // no vectors to be orthogonalized
  {
    return;
  }
  if (R1==NULL)
  {
    // return "bad cast or null pointer"
    *ierr=PHIST_BAD_CAST;
    return;
  }

  MT normW0[k];
  MT normW1[k];

  if (m==0)
  {
    SUBR(mvec_QR)(W,R1,ierr);
    if (*ierr==PHIST_NOT_IMPLEMENTED)
    {
      PHIST_CHK_NEG_IERR(SUBR(mvec_QR_fallback)(W,R1,normW1,ierr),*ierr);
      // indicate that R1 is invalid.
      // note: information about the rank of V is lost!
      *ierr=-7;
    }
    else
    {
      PHIST_CHK_NEG_IERR(*ierr,*ierr);
    }
    return;
  }

  if (R2==NULL)
  {
    // return "bad cast or null pointer"
    *ierr=PHIST_BAD_CAST;
    return;
  }

      PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(sdMat_create)(&R1p,k,k,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R1pp,k,k,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R2p,m,k,comm,ierr),*ierr);


  // check that all array dimensions are correct
  PHIST_CHK_IERR(*ierr=(int)(V==NULL || W==NULL || R1==NULL || R2==NULL),*ierr);

  int nrR1,ncR1,nrR2,ncR2;
  lidx_t n, tmp;

  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(W,&tmp,ierr),*ierr);

  PHIST_CHK_IERR(*ierr=((n==tmp)?0:-1),*ierr);

  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R1,&nrR1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R1,&ncR1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R2,&nrR2,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R2,&ncR2,ierr),*ierr);

  // Q (and W) must match R1, R1 must be square
  PHIST_CHK_IERR(*ierr=((k==nrR1)?0:-1),*ierr);
  PHIST_CHK_IERR((*ierr=(k==ncR1)?0:-1),*ierr);
  // V must match R2, and V*R2 must match W
  PHIST_CHK_IERR((*ierr=(m==nrR2)?0:-1),*ierr);
  PHIST_CHK_IERR((*ierr=(k==ncR2)?0:-1),*ierr);

#ifdef TESTING
  PHIST_DEB("orthog: V is %dx%d,  W is %dx%d\n"
                        "       R1 is %dx%d, R2 is %dx%d\n",n,m,n,k,nrR1,ncR1,nrR2,ncR2);
#endif

  // compute the norms of the columns in W (for checking the result later)
  
  // determine original norms of W vectors and breakdown tolerance
  PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW0,ierr),*ierr);
  breakdown = normW0[0];
  PHIST_DEB("orthog: normW0 is");
  for (int i=0;i<k;i++)
  {
    breakdown=std::min(normW0[i], breakdown);
    PHIST_DEB(" %e", normW0[i]);
  }
  breakdown*=mt::eps()*1000;
  PHIST_DEB("\n");
  
  // orthogonalize against V (first CGS sweep)
#ifdef TESTING
{
  _MT_ normV[m];
  PHIST_CHK_IERR(SUBR(mvec_norm2)(V,normV,ierr),*ierr);
  PHIST_DEB("orthog: normV is");
  for(int i = 0; i < m; i++)
  {
    PHIST_DEB(" %e", normV[i]);
  }
  PHIST_DEB("\n");
}
#endif



  //R2=V'*W;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2,ierr),*ierr);
#ifdef TESTING
  PHIST_SOUT(PHIST_INFO,"R2:\n");
  PHIST_CHK_IERR(SUBR(sdMat_print)(R2, ierr), *ierr);
#endif
  //W=W-V*R2;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,R2,st::one(),W,ierr),*ierr);
  // norms after first CGS sweep. This could be done cheaper by using R1 computed below,
  // since ||W_j||_2 = ||R1_j||_2 (TODO - that's how it's done in Belos).
  PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW1,ierr),*ierr);
  PHIST_DEB("orthog: normW1[0] is %e\n", normW1[0]);

  // special case which cannot be detected in mvec_QR if W in span(V) (because mvec_QR uses a relative tolerance):
  MT maxNormW1 = 0;
  for(int i = 0; i < k; i++)
    maxNormW1 = std::max(normW1[i],maxNormW1);
  if(maxNormW1 < breakdown)
  {
    PHIST_SOUT(PHIST_INFO,"breakdown in phist_orthog: W in span(V), filling in random vectors!\n");
    PHIST_CHK_IERR(SUBR(mvec_put_value)(W,st::zero(),ierr),*ierr);
  }
  // orthogonalize W after the first GS pass. This gives us the rank  
  // of W (ierr>0 indicates the dimension of the null space of W, i.e.
  // the rank of W is k-ierr). If W is not full rank, it is augmented 
  // by random vectors. These still have to be orthogonalized against 
  // V, which we do next. The projection coefficients for this part   
  // are thrown away as the randomized vectors are not really related 
  //to W anyway.
  SUBR(mvec_QR)(W,R1,ierr);
  if (*ierr!=PHIST_NOT_IMPLEMENTED)
  {
    PHIST_CHK_NEG_IERR(*ierr,*ierr);
    rankW=k-*ierr;
    // set normW1 to diag(R1)
    _ST_ *R1_raw = NULL;
    lidx_t ldaR1;
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(R1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R1,&R1_raw,&ldaR1,ierr),*ierr);
    for(int i = 0; i < k; i++)
    {
      normW1[i] = std::min(normW1[i],st::abs(R1_raw[i+ldaR1*i]));
    }
  }
  else
  {
    // recomputes normR1 along the way, so we might as well use the memory
    // allocated for that beforehand (even though the result is already there)
    PHIST_CHK_NEG_IERR(SUBR(mvec_QR_fallback)(W,R1,normW1,ierr),*ierr);
    rankW=k-*ierr;
    useSVQB=true;
  }
#ifdef TESTING
  PHIST_SOUT(PHIST_VERBOSE,"rankW is %d, R1:\n", rankW);
  PHIST_CHK_IERR(SUBR(sdMat_print)(R1, ierr), *ierr);
  if (useSVQB)
  {
    PHIST_SOUT(PHIST_VERBOSE,"orthog uses fallback SVQB kernel\n");
  }
#endif
  *ierr = 0;


  if (rankW < k )
  {
    int random_iter = 0;
    do
    {
      // terminate even if random vectors are not "random" enough, e.g. random number generator is broken
      if( random_iter++ > 10 )
        {
        PHIST_SOUT(PHIST_ERROR,"could not create random orthogonal vectors, possibly the random vector generator is broken!\n");
        *ierr = -8;
        return;
        }
      st::sdMat_t *Rrnd;
      PHIST_SOUT(PHIST_INFO,"Matrix W does not have full rank (%d cols, rank=%d)\n",k,rankW);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&Rrnd,m,k,comm,ierr),*ierr);
      //R2=V'*W;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),Rrnd,ierr),*ierr);
      //W=W-V*R2;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,Rrnd,st::one(),W,ierr),*ierr);
      // throw away projection coefficients.
      PHIST_CHK_IERR(SUBR(sdMat_delete)(Rrnd,ierr),*ierr);

      // reorthogonalize result, filling in new random vectors if these were in span(V)
      if (useSVQB)
      {
        MT dum[k];
        PHIST_CHK_NEG_IERR(SUBR(mvec_QR_fallback)(W,R1,dum,ierr),*ierr);
      }
      else
      {
        PHIST_CHK_NEG_IERR(SUBR(mvec_QR)(W,R1,ierr),*ierr);
      }
    } while (*ierr > 0);

    // we must fill the appropriate columns of R1 with zeros (where random values were used)
    TYPE(sdMat_ptr) R1_r = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(R1,&R1_r,0,k-1,rankW,k-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R1_r,st::zero(),ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_r,ierr),*ierr);
  }
    
  int step=1;
  while (true)
  {
    MT maxRed=1.0;
    for (int j=0;j<k;j++)
    {
      maxRed=std::min(maxRed,normW1[j]/normW0[j]);
      normW0[j]=1.; // after mvec_QR the norm is 1!
    }
    PHIST_SOUT(PHIST_VERBOSE,"reduction in norm, GS step %d: %4.2f\n",step,maxRed);
    if (maxRed>0.85)
      {
      stopGS=true;
      PHIST_SOUT(PHIST_VERBOSE,"stopping Gram-Schmidt\n");
      break;
      }
    if (step>=numSweeps)
      {
      PHIST_SOUT(PHIST_VERBOSE,"stopping Gram-Schmidt because %d steps have been performed\n", 
                numSweeps);
      }
    step++;
    //R2p=V'*W;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2p,ierr),*ierr);

    //W=W-V*R2';
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,R2p,st::one(),W,ierr),*ierr);

    // we must not modify columns in R2 corresponding to random vectors!
    if( rankW < k )
    {
      TYPE(sdMat_ptr) R2p_rand = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(R2p,&R2p_rand,0,m-1,rankW,k-1,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R2p_rand,st::zero(),ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p_rand,ierr),*ierr);
    }

    //R2=R2+R2';
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R2p,st::one(),R2,ierr),*ierr);

    // again, the norm can be computed from the coefficients in R2p (TODO)
    PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW1,ierr),*ierr);    

    // orthogonalize W. The situation where the rank of W becomes
    // smaller than k here is (I think) very unlikely because we 
    // added random vectors before if rank(W) was not full. If   
    // it happens, we return with error code -10.                
    if (useSVQB)
    {
      PHIST_CHK_NEG_IERR(SUBR(mvec_QR_fallback)(W,R1p,normW1,ierr),*ierr);
      // we do not update R1, it will just be flagged invalid at the end
      // of the subroutine by setting *ierr=-7
    }
    else
    {
      PHIST_CHK_NEG_IERR(SUBR(mvec_QR)(W,R1p,ierr),*ierr);

      if (*ierr>0)
      {
        PHIST_SOUT(PHIST_ERROR,"Unexpected rank deficiency in orthog routine\n(file %s, line %d)\n",
                __FILE__,__LINE__);
        *ierr=-10;
        return;
      }

      _ST_ *R1p_raw = NULL;
      lidx_t ldaR1p;
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(R1p,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R1p,&R1p_raw,&ldaR1p,ierr),*ierr);
      for(int i = 0; i < k; i++)
        normW1[i] = std::min(normW1[i],st::abs(R1p_raw[i+ldaR1p*i]));

      //R1pp=R1
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1pp,ierr),*ierr);
      //R1=R1p*R1pp;
      PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R1p,R1pp,st::zero(),R1,ierr),*ierr);
      // keep zero entries in R1 for rank(W-V*V'*W) < k
      if( rankW < k )
      {
        // we must fill the appropriate columns of R1 with zeros (where random values were used)
        TYPE(mvec_ptr) R1_r = NULL;
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(R1,&R1_r,0,k-1,rankW,k-1,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(sdMat_put_value)(R1_r,st::zero(),ierr),*ierr);
        PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_r,ierr),*ierr);
      }
    }//mvec_QR available?
  }//while

  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1p,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1pp,ierr),*ierr);

  // return size of randomly filled null space if W-V*R2 not full rank
  if (rankW < k) *ierr = k-rankW;

  // if in the last CGS sweep a column decreased too much in norm, 
  // return an error (we could return a warning, but since we used 
  // all positive numbers for the case of rank deficient W already,
  // the caller will have to deal with the -5 case specially.      
  if (stopGS==false)
  {
    // we need to check again, for instance, the 2nd pass may have given
    // the desired small reduction but the check would have been performed
    // at the beginning of pass 3, which was not done.
    for (int i=0;i<k;i++)
    {
      if (normW1[i]<0.7*normW0[i])
      {
        // more orthogonalization steps may be required.
        *ierr=-9;
      }
    }
  }//stopGS?
  
  if (useSVQB && ierr>=0) *ierr=-7;
  return;
}

void SUBR(mvec_QR_fallback)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, _MT_* nrmsV, int *ierr)
{
      PHIST_ENTER_FCN(__FUNCTION__);
      int k;
      PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&k,ierr),*ierr);
      _MT_ dummy[k];
      _MT_* nrmsV_ptr=nrmsV;
      int dim0=k;
      int return_value;
      // try a couple of times to get a full rank
      // orthogonal matrix, otherwise just give up
      for (int i=0;i<3;i++)
      {
        // use fallback kernel: SVQB
        PHIST_CHK_NEG_IERR(SUBR(svqb)(V,B,nrmsV,ierr),*ierr);
        // next sweep do not overwrite nrmsV
        nrmsV_ptr=dummy;
        dim0=*ierr;
        if (i==0) return_value=dim0;
        if (dim0==0) break;
        TYPE(mvec_ptr) V0=NULL;
        PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&V0,0,dim0-1,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_random)(V0,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_delete)(V0,ierr),*ierr);
      }
      *ierr=return_value;// dimension of nullspace found infirst sweep
      return;
}
