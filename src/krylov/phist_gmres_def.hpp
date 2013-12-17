//! small helper function for GMRES to compute
// a Givens rotation such that
// | c s ||f| |r|
// | _   || |=| |
// |-s c ||g| |0|
// This function ses the convention for the BLAS,
// cf. for instance Bindel et al. 2001). The
// BLAS routine is called XROTG, however, it does
// not do what I want in the complex case (TODO is there
// another routine that gives complex coefficients??)
void SUBR(rotg)(_ST_ f, _ST_ g, _ST_& cs, _ST_& sn, _ST_& r)
{
#include "phist_std_typedefs.hpp"
  MT af=st::abs(f);
  MT ag=st::abs(g);
  MT d=mt::sqrt(af*af + ag*ag);
  cs=st::zero();
  sn=st::zero();
  r=st::zero();

  if (ag==mt::zero() ) // includes case f=0 and g=0
  {
    cs = st::one();
    sn = st::zero();
    r  = f;
  }
  else if (af == mt::zero())
  {
    cs=st::zero();
    sn=st::conj(g)/ag;
    r=ag;
  }
  else
  {
    cs=af/d;
    sn=(f/af)*(st::conj(g)/d);
    r=(f/af)*d;
  }
  return;
}

// create new state object
void SUBR(gmresState_create)(TYPE(gmresState_ptr)* state, const_map_ptr_t map, 
        int maxBas,int* ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  TYPE(gmresState)* S = new TYPE(gmresState);
  *state = S;
  // options not set
  S->id=-1;
  S->tol=-mt::one(); 
  S->maxIters=-1;
  S->maxBas=maxBas;
  S->ierr=-2;

  S->maxBasAllocated_=maxBas;

  PHIST_CHK_IERR(SUBR(mvec_create)(&S->X0_, map,1, ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&S->V_, map,maxBas, ierr),*ierr);
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S->H_, maxBas+1, maxBas, comm,ierr),*ierr);
  S->cs_ = new ST[maxBas];
  S->sn_ = new ST[maxBas];
  S->rs_ = new ST[maxBas];
  S->curDimV_=0;
  S->curIter_=0;
  S->normB_=-mt::one(); // not initialized
}

//! delete gmresState object
void SUBR(gmresState_delete)(TYPE(gmresState)* S, int* ierr)
{
  PHIST_CHK_IERR(SUBR(mvec_delete)(S->V_,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(S->H_,ierr),*ierr);
  delete [] S->cs_;
  delete [] S->sn_;
  delete [] S->rs_;
  delete S;
}

// reset gmres state.
void SUBR(gmresState_reset)(TYPE(gmresState_ptr) S, TYPE(const_mvec_ptr) b,
        TYPE(const_mvec_ptr) x0,int *ierr)
{
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  if (b==NULL && S->B_==NULL)
    {
    PHIST_OUT(PHIST_ERROR,"on the first call to gmresState_reset you *must* provide the RHS vector");
    *ierr=-1;
    return;
    }
  else if (b!=S->B_ && b!=NULL)
    {
    // new rhs -> need to recompute ||B||
    S->B_=b;
    S->normB_=-mt::one();
    S->curIter_=0;
    }
  S->curDimV_=0;
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),x0,st::zero(),S->X0_,ierr),*ierr);
  return;
}

void SUBR(gmresState_updateSol)(TYPE(gmresState_ptr) S_array[], TYPE(mvec_ptr) x, int* ierr)
{
#include "phist_std_typedefs.hpp"

  int numSys;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&numSys,ierr),*ierr);
  const_comm_ptr_t comm=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_comm)(x,&comm,ierr),*ierr);

  // TODO - omp parallel for
  for (int i=0;i<numSys;i++)
  {
    // compute y by solving the triangular system
    TYPE(const_gmresState_ptr) S = S_array[i];
    TYPE(sdMat_ptr) y=NULL;
    ST* H_raw=NULL,*y_raw=NULL;
    lidx_t ldH,ldy;
    
    PHIST_CHK_IERR(SUBR(sdMat_create)(&y,S->curDimV_,numSys,comm,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(y,&y_raw,&ldy,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(S->H_,&H_raw,&ldH,ierr),*ierr);

    // y = H\rs, H upper triangular
    const char* uplo="U";
    const char* trans="N";
    const char* diag="N";
    int nrhs=1;
    PHIST_CHK_IERR(PREFIX(TRTRS)(uplo,trans,diag,&S->curDimV_,&nrhs,
                                        (st::blas_scalar_t*)H_raw,&ldH,
                                        (st::blas_scalar_t*)y_raw, &ldy, ierr),*ierr);

    // X = X + M\V*y. TODO: with right preconditioning, split this into two parts and apply
    // preconditioner to all systems simultaneously outside the loop (need tmp vector)
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),S->V_,y,st::one(),x,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(y,ierr),*ierr);
  }
  return;
}

// implementation of gmres on several systems simultaneously
void SUBR(gmresState_iterate)(TYPE(const_op_ptr) Op,
        TYPE(gmresState_ptr) S_array[], int numSys,
        int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // map defining how to create new vector objects
  const_map_ptr_t map=NULL;
  // multi-vectors to work with several of the systems simultaneously
  TYPE(mvec_ptr) R=NULL, V=NULL, W=NULL;
  // views into vector blocks of system i
  TYPE(mvec_ptr) Vprev=NULL, Vj=NULL;
  // for the orthog routine (again, a view into the H objects in the states
  TYPE(sdMat_ptr) R1=NULL,R2=NULL;
  
  // check if the given state objects have computed the norm of B,
  // and do so if not.
  for (int i=0;i<numSys;i++)
  {
    S_array[i]->ierr=1; // not converged yet
    if (S_array[i]->normB_<mt::zero())
      {
        SUBR(mvec_norm2)(S_array[i]->B_,&S_array[i]->normB_,&S_array[i]->ierr);
        PHIST_CHK_IERR(*ierr=S_array[i]->ierr,*ierr);
      }
  }
  
  // figure out how new vectors are built
  PHIST_CHK_IERR(SUBR(mvec_get_map)(S_array[0]->X0_,&map,ierr),*ierr);

//TODO - we could use "scattered views" here to avoid copying the columns into the GMRES 
//       states

  // residual vectors
  PHIST_CHK_IERR(SUBR(mvec_create)(&R,map,numSys,ierr),*ierr);
  // vectors V(:,j) for each of the systems
  PHIST_CHK_IERR(SUBR(mvec_create)(&V,map,numSys,ierr),*ierr);
  // vectors W=A*V(:,j) for each of the systems
  PHIST_CHK_IERR(SUBR(mvec_create)(&W,map,numSys,ierr),*ierr);

  bool anyConverged=false;
  bool anyFailed=false;

//TODO - here we must compute the initial residual r0=A*x0-b, rs0=||r0||, v0=r0/rs0 after
//       a call to reset(), but without interfering with the rest of the systems...

  // we return as soon as one system converges or reaches its
  // maximum permitted number of iterations. The decision about what to do
  // next is then left to the caller.

  // Arnoldi - build orthogonal basis V and upper Hessenberg matrix H.
  // jmax is chosen so that the next system reaches maxBas.
  int jmax = S_array[0]->maxBas-S_array[0]->curDimV_;
  for (int i=1;i<numSys;i++)
  {
    jmax = std::min(jmax, S_array[i]->maxBas-S_array[i]->curDimV_);
  }
  for (int j_dum=0; j_dum<jmax; j_dum++)
  {
    //W=A*(M\V(:,j));
    PHIST_CHK_IERR(Op->apply(st::one(),Op->A, V, st::zero(), W, ierr), *ierr);
      
    // Arnoldi update. TODO - we could save some messages by
    // clustering the communication in orthog
    // TODO - maybe we could parallelize this loop by an OpenMP section
    //        so that the small stuff can be done in parallel too?
    for (int i=0;i<numSys;i++)
    {
      TYPE(gmresState_ptr) S = S_array[i];
      int j=S->curDimV_;
      
      if (j<=0)
        {
        PHIST_OUT(PHIST_ERROR,"gmres state not initialized correctly");
        S->ierr=-1;
        *ierr=-1;
        }
      if (j>=S->maxBas-1)
        {
        PHIST_OUT(PHIST_ERROR,"gmres state not initialized/reset correctly");
        S->ierr=2;
        *ierr=2;
        }
        
      //TODO - blocking of vectors in orthog to avoid communication??
        
      // view the existing basis vectors as Vprev
      PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Vprev,0,j-1,ierr),*ierr);
      // view the next V vector as Vj
      PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Vj,j,j,ierr),*ierr);
      // copy vector A*v into location Vj
      PHIST_CHK_IERR(SUBR(mvec_get_block)(V,Vj,i,i,ierr),*ierr);
      // view H(j,j) as R1
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_,&R1,j,j,j,j,ierr),*ierr);
      // view H(1:j,j) as R2
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_,&R2,0,j-1,j,j,ierr),*ierr);
      //orthogonalize
      PHIST_CHK_IERR(SUBR(orthog)(Vprev,Vj,R1,R2,2,ierr),*ierr);
        
      //    % update QR factorization of H

      // note: it is OK to work on raw data of serial dense matrices because 
      // they are typically not modified on an accelerator, however, we must 
      // somehow implement a check of this in the kernel interfaces (TODO)   
      // note also that we will access Hcol[j], which is strictly speaking   
      // R1, but by construction they should be aligned in memory. 
      ST *Hcol=NULL; 
      lidx_t ldH; 
      PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R2,&Hcol,&ldH,ierr),*ierr); 

      ST htmp;

      //    % apply previous (j-1) transformations to column j
      for (int jj=0;jj<j-1;jj++)
      {
        htmp = S->cs_[jj]*Hcol[jj] + 
                    S->sn_[jj]*Hcol[jj+1]; // H(j-1,j) in the last step, which is R1 above
        Hcol[jj+1] = -st::conj(S->sn_[jj])*Hcol[jj] + S->cs_[jj]*Hcol[jj+1];
        Hcol[jj]   = htmp;
      }

      // new Givens rotation for eliminating H(j+1,j)
      SUBR(rotg)(Hcol[j-1],Hcol[j],S->cs_[j-1],S->sn_[j-1],htmp);
      // eliminate H(j,j-1)
      Hcol[j-1] = htmp;
      Hcol[j]=st::zero();
      // apply to RHS
      htmp=S->cs_[j-1]*S->rs_[j-1];
      S->rs_[j] = -st::conj(S->sn_[j-1])*S->rs_[j-1];
      S->rs_[j-1]=htmp;

      if (PHIST_OUTLEV>=PHIST_DEBUG)
      {
        PHIST_DEB("transformed H");
        PHIST_CHK_IERR(SUBR(sdMat_print)(S->H_,ierr),*ierr);
        //disp(rs(1:j+1,:))
      }
      S->normR_=st::abs(S->rs_[j+1]);
      MT relres=S->normR_/S->normB_;
      
      if (relres<S->tol)
      {
        S->ierr=0; // converged
        anyConverged=true;
      }
    }// for all systems i, update H and rs
    if (anyConverged) break;
  }// for-loop (Arnoldi iterations)
 
 
  if (anyConverged)
  {
    *ierr=0;
  }
      
  if (anyFailed)
  {
    *ierr=1;
  }
return;
}
