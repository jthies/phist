#ifdef HAVE_TRILINOS_ORTHO_MANAGER
// provide trili_orthog as fallback if B!=NULL
#include "trili_orthog_def.hpp"
#endif
// helper routine in case the kernel lib does not provide mvec_QR, which
// is quite an advanced kernel to ask for.
void SUBR(mvec_QB)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, _MT_* nrmsV, int *iflag);

//! orthogonalize an mvec against an already orthogonal one.
void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_op_ptr) B,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps,
                     int* rankVW,
                     int* iflag)
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

  if (B!=NULL)
  {
#ifdef HAVE_TRILINOS_ORTHO_MANAGER
  static bool first_call=true;
  if (first_call)
  {
    PHIST_SOUT(PHIST_WARNING,"using trilinos ortho manager for B-orthogonalization\n");
    first_call=false;
  }
  
  SUBR(trili_orthog)(V,W,B,R1,R2,numSweeps,rankVW,iflag);
  return;    
#else
    PHIST_SOUT(PHIST_ERROR,"case B!=I not implemented (file %s, line %d)\n",__FILE__,__LINE__);
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
#endif
  }

  bool useSVQB=false;
  *iflag=0;
  *rankVW=-1;
  
  if (V==NULL)
  {
    m=0;
  }
  else
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
  }
  if (W==NULL)
  {
    k=0;
  }
  else
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
  }
  if (k==0) // no vectors to be orthogonalized
  {
    *rankVW=m;
    return;
  }
  if (R1==NULL)
  {
    // return "bad cast or null pointer"
    *iflag=PHIST_BAD_CAST;
    return;
  }

  MT normW0[k];
  MT normW1[k];

  if (m==0)
  {
    SUBR(mvec_QR)(W,R1,iflag);
    if (*iflag==PHIST_NOT_IMPLEMENTED)
    {
      PHIST_CHK_NEG_IERR(SUBR(mvec_QB)(W,R1,normW1,iflag),*iflag);
      *rankVW=k-*iflag;
      // indicate that R1 is invalid.
      *iflag=+2;
    }
    else
    {
      PHIST_CHK_NEG_IERR((void)*iflag,*iflag);
      *rankVW=k-*iflag;
      *iflag=std::min(*iflag,+1); // iflag=+1 means: [V,W] is rank deficient
    }
    return;
  }

  if (R2==NULL)
  {
    // return "bad cast or null pointer"
    *iflag=PHIST_BAD_CAST;
    return;
  }

      PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_create)(&R1p,k,k,comm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R1pp,k,k,comm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&R2p,m,k,comm,iflag),*iflag);


  // check that all array dimensions are correct
  PHIST_CHK_IERR(*iflag=(int)(V==NULL || W==NULL || R1==NULL || R2==NULL),*iflag);

  int nrR1,ncR1,nrR2,ncR2;
  lidx_t n, tmp;

  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(W,&tmp,iflag),*iflag);

  PHIST_CHK_IERR(*iflag=((n==tmp)?0:-1),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R1,&nrR1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R1,&ncR1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R2,&nrR2,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R2,&ncR2,iflag),*iflag);

  // Q (and W) must match R1, R1 must be square
  PHIST_CHK_IERR(*iflag=((k==nrR1)?0:-1),*iflag);
  PHIST_CHK_IERR((*iflag=(k==ncR1)?0:-1),*iflag);
  // V must match R2, and V*R2 must match W
  PHIST_CHK_IERR((*iflag=(m==nrR2)?0:-1),*iflag);
  PHIST_CHK_IERR((*iflag=(k==ncR2)?0:-1),*iflag);

#ifdef TESTING
  PHIST_DEB("orthog: V is %dx%d,  W is %dx%d\n"
                        "       R1 is %dx%d, R2 is %dx%d\n",n,m,n,k,nrR1,ncR1,nrR2,ncR2);
#endif

  // compute the norms of the columns in W (for checking the result later)
  
  // determine original norms of W vectors and breakdown tolerance
  PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW0,iflag),*iflag);
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
  PHIST_CHK_IERR(SUBR(mvec_norm2)(V,normV,iflag),*iflag);
  PHIST_DEB("orthog: normV is");
  for(int i = 0; i < m; i++)
  {
    PHIST_DEB(" %e", normV[i]);
  }
  PHIST_DEB("\n");
}
#endif



  //R2=V'*W;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2,iflag),*iflag);
#ifdef TESTING
  //PHIST_SOUT(PHIST_INFO,"R2:\n");
  //PHIST_CHK_IERR(SUBR(sdMat_print)(R2, iflag), *iflag);
#endif
  //W=W-V*R2;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,R2,st::one(),W,iflag),*iflag);
  // norms after first CGS sweep. This could be done cheaper by using R1 computed below,
  // since ||W_j||_2 = ||R1_j||_2. Since this is assuming exact arithmetic, we leave it
  // for now and check for a (near) breakdown before the normalization step. The 
  // normalization (mvec_QR or SVQB) will fill in random numbers in that case.
  PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW1,iflag),*iflag);
  PHIST_DEB("orthog: normW1[0] is %e\n", normW1[0]);

  // special case which cannot be detected in mvec_QR if W in span(V) (because mvec_QR uses a relative tolerance):
  MT maxNormW1 = 0;
  for(int i = 0; i < k; i++)
    maxNormW1 = std::max(normW1[i],maxNormW1);
  if(maxNormW1 < breakdown)
  {
    PHIST_SOUT(PHIST_INFO,"breakdown in phist_orthog: W in span(V), filling in random vectors!\n");
    PHIST_CHK_IERR(SUBR(mvec_put_value)(W,st::zero(),iflag),*iflag);
  }
  // orthogonalize W after the first GS pass. This gives us the rank  
  // of W (iflag>0 indicates the dimension of the null space of W, i.e.
  // the rank of W is k-iflag). If W is not full rank, it is augmented 
  // by random vectors. These still have to be orthogonalized against 
  // V, which we do next. The projection coefficients for this part   
  // are thrown away as the randomized vectors are not really related 
  //to W anyway.
  SUBR(mvec_QR)(W,R1,iflag);
  if (*iflag!=PHIST_NOT_IMPLEMENTED)
  {
    PHIST_CHK_NEG_IERR((void)*iflag,*iflag);
    rankW=k-*iflag;
    // set normW1 to diag(R1)
    _ST_ *R1_raw = NULL;
    lidx_t ldaR1;
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(R1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R1,&R1_raw,&ldaR1,iflag),*iflag);
    for(int i = 0; i < k; i++)
    {
      normW1[i] = std::min(normW1[i],st::abs(R1_raw[i+ldaR1*i]));
    }
  }
  else
  {
    // recomputes normR1 along the way, so we might as well use the memory
    // allocated for that beforehand (even though the result is already there)
    PHIST_CHK_NEG_IERR(SUBR(mvec_QB)(W,R1,normW1,iflag),*iflag);
    rankW=k-*iflag;
    useSVQB=true;
  }
#ifdef TESTING
  //PHIST_SOUT(PHIST_VERBOSE,"rankW is %d, R1:\n", rankW);
  //PHIST_CHK_IERR(SUBR(sdMat_print)(R1, iflag), *iflag);
  if (useSVQB)
  {
    PHIST_SOUT(PHIST_VERBOSE,"orthog uses fallback SVQB kernel\n");
  }
#endif
  *iflag = 0;


  if (rankW < k )
  {
    int random_iter = 0;
    do
    {
      // terminate even if random vectors are not "random" enough, e.g. random number generator is broken
      if( random_iter++ > 10 )
        {
        PHIST_SOUT(PHIST_ERROR,"could not create random orthogonal vectors, possibly the random vector generator is broken!\n");
        *iflag = -8;
        return;
        }
      st::sdMat_t *Rrnd;
      PHIST_SOUT(PHIST_INFO,"Matrix W does not have full rank (%d cols, rank=%d)\n",k,rankW);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&Rrnd,m,k,comm,iflag),*iflag);
      //R2=V'*W;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),Rrnd,iflag),*iflag);
      //W=W-V*R2;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,Rrnd,st::one(),W,iflag),*iflag);
      // throw away projection coefficients.
      PHIST_CHK_IERR(SUBR(sdMat_delete)(Rrnd,iflag),*iflag);

      // reorthogonalize result, filling in new random vectors if these were in span(V)
      if (useSVQB)
      {
        MT dum[k];
        PHIST_CHK_NEG_IERR(SUBR(mvec_QB)(W,R1,dum,iflag),*iflag);
      }
      else
      {
        PHIST_CHK_NEG_IERR(SUBR(mvec_QR)(W,R1,iflag),*iflag);
      }
    } while (*iflag > 0);

    // we must fill the appropriate columns of R1 with zeros (where random values were used)
    TYPE(sdMat_ptr) R1_r = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(R1,&R1_r,0,k-1,rankW,k-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R1_r,st::zero(),iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_r,iflag),*iflag);
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
      break;
    }
    step++;
    //R2p=V'*W;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),R2p,iflag),*iflag);

    //W=W-V*R2';
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),V,R2p,st::one(),W,iflag),*iflag);

    // we must not modify columns in R2 corresponding to random vectors!
    if( rankW < k )
    {
      TYPE(sdMat_ptr) R2p_rand = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(R2p,&R2p_rand,0,m-1,rankW,k-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R2p_rand,st::zero(),iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p_rand,iflag),*iflag);
    }

    //R2=R2+R2';
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R2p,st::one(),R2,iflag),*iflag);

    // again, the norm can be computed from the coefficients in R2p (TODO)
    PHIST_CHK_IERR(SUBR(mvec_norm2)(W,normW1,iflag),*iflag);    

    // orthogonalize W. The situation where the rank of W becomes
    // smaller than k here is (I think) very unlikely because we 
    // added random vectors before if rank(W) was not full. If   
    // it happens, we return with error code -10.                
    if (useSVQB)
    {
      PHIST_CHK_NEG_IERR(SUBR(mvec_QB)(W,R1p,normW1,iflag),*iflag);
      // we do not update R1, it will just be flagged invalid at the end
      // of the subroutine by setting *iflag=+2
    }
    else
    {
      PHIST_CHK_NEG_IERR(SUBR(mvec_QR)(W,R1p,iflag),*iflag);

      if (*iflag>0)
      {
        PHIST_SOUT(PHIST_ERROR,"Unexpected rank deficiency in orthog routine\n(file %s, line %d)\n",
                __FILE__,__LINE__);
        *iflag=-10;
        return;
      }

      _ST_ *R1p_raw = NULL;
      lidx_t ldaR1p;
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(R1p,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R1p,&R1p_raw,&ldaR1p,iflag),*iflag);
      for(int i = 0; i < k; i++)
        normW1[i] = std::min(normW1[i],st::abs(R1p_raw[i+ldaR1p*i]));

      //R1pp=R1
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1pp,iflag),*iflag);
      //R1=R1p*R1pp;
      PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R1p,R1pp,st::zero(),R1,iflag),*iflag);
      // keep zero entries in R1 for rank(W-V*V'*W) < k
      if( rankW < k )
      {
        // we must fill the appropriate columns of R1 with zeros (where random values were used)
        TYPE(mvec_ptr) R1_r = NULL;
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(R1,&R1_r,0,k-1,rankW,k-1,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_put_value)(R1_r,st::zero(),iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_r,iflag),*iflag);
      }
    }//mvec_QR available?
  }//while

  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1p,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1pp,iflag),*iflag);

  if (rankW < k)
  { 
    *iflag = +1;
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
    for (int i=0;i<k;i++)
    {
      if (normW1[i]<0.7*normW0[i])
      {
        // more orthogonalization steps may be required.
        *iflag=-9;
      }
    }
  }//stopGS?

  *rankVW=m+rankW;
  
  if (useSVQB) *iflag=+2;
  return;
}

void SUBR(mvec_QB)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, _MT_* nrmsV, int *iflag)
{
      PHIST_ENTER_FCN(__FUNCTION__);
      int k;
      PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&k,iflag),*iflag);
      _MT_ dummy[k];
      _MT_* nrmsV_ptr=nrmsV;
      int dim0=k;
      int return_value;
      // try a couple of times to get a full rank
      // orthogonal matrix, otherwise just give up.
      // We do at least 2 sweeps for robustness' sake.
      for (int i=0;i<3;i++)
      {
        // use fallback kernel: SVQB
        PHIST_CHK_NEG_IERR(SUBR(svqb)(V,B,nrmsV_ptr,iflag),*iflag);
        // next sweep do not overwrite nrmsV
        nrmsV_ptr=dummy;
        dim0=*iflag;
        if (i==0) return_value=dim0;
        else if (dim0==0) break;
        TYPE(mvec_ptr) V0=NULL;
        PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&V0,0,dim0-1,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_random)(V0,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_delete)(V0,iflag),*iflag);
      }
      *iflag=return_value;// dimension of nullspace found infirst sweep
      return;
}
