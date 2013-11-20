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
  ST d=st::sqrt(std::abs(f).^2 + abs(g).^2);
  c=st::zero();
  s=c; r=c;

  if (g==st::zero() ) // includes case f=0 and g=0
    {
    cs = st::one();
    sn = zt::zero();
    r  = f;
    }
  else if (f == st::zero())
    {
    c=st::zero();
    s=st::conj(g)/std::abs(g);
    r=std::abs(g);
    }
  else
    {
    c=std::abs(f)/d;
    s=(f/std::abs(f))*(st::conj(g)/d);
    r=(f/std::abs(f))*d;
    }
  }

// create new state object
void SUBR(gmresState_create)(TYPE(gmresState_ptr)* state, const_map_ptr map, 
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
  PHIST_CHK_IERR(SUBR(mvec_create)(&S->V_, maxBas, ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&S->H_, maxBas+1, maxBas, ierr),*ierr);
  S->cs_ = new ST[maxBas];
  S->sn_ = new ST[maxBas];
  S->rs_ = new ST[maxBas];
  S->curDim_=0;
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
SUBR(gmresState_reset)(TYPE(gmresState_ptr) S, TYPE(const_mvec_ptr) b,
        TYPE(const_mvec_ptr) x0,int *ierr)
  {
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  MT rs0;
  TYPE(mvec_ptr) v0 = NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V_,&v0,0,0,ierr),*ierr);
  if (b!=NULL)
    {
    PHIST_CHK_IERR(SUBR(mvec_norm2)(b,&(S->normB_),ierr),*ierr);
    }
  else if (S->normB_==-mt::one())
    {
    PHISt_OUT(PHIST_ERROR,"on the first call to gmresState_reset you *must* provide the RHS vector");
    *ierr=-1;
    return;
    }
  if (x0!=NULL)
    {
    // copy starting vector into V
    PHIST_CHK_IERR(SUBR(mvec_set_block)(S->V_,x0,0,0,ierr),*ierr);
    }
  else
    {
    // generate random initial vector
    PHIST_CHK_IERR(SUBR(mvec_random)(v0,ierr),*ierr);
    }
  //normalize
  PHIST_CHK_IERR(SUBR(mvec_normalize)(v0,&rs0,ierr),*ierr);
  S->rs_[0]=(ST)rs0;
  S->curIter_=1;
  S->curDim_=1;
  // delete view of first column
  PHIST_CHK_IERR(mvec_delete)(v0,ierr),*ierr);
  return;
  }

// implementation of gmres on several systems simultaneously
void SUBR(gmres)(TYPE(const_op_ptr) Op,
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B,
        TYPE(gmresState)* S_array,
        int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // map defining how to create new vector objects
  const_map_t* map=NULL;
  // multi-vectors to work with several of the systems simultaneously
  TYPE(mvec_ptr) R=NULL, V=NULL, W=NULL;
  // views into vector blocks of system i
  TYPE(mvec_ptr) Vprev=NULL, Vj=NULL;
  // for the orthog routine (again, a view into the H objects in the states
  TYPE(sdMat_ptr) R1=NULL,R2=NULL;
  
  // number of systems we're iterating on. This value may decrease
  // as systems converge, until at least 'stopIf' of them are
  // converged and we exit the subroutine.
  int numSys;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&numSys,ierr),*ierr);
  
  // check if the given state objects have computed the norm of B,
  // and do so if not.
  for (int i=0;i<numSys;i++)
    {
    S_array[i]->ierr_=1; // not converged yet
    }
  
  MT beta[numSys];
  
  // figure out how new vectors are built
  PHIST_CHK_IERR(SUBR(mvec_get_map)(X,&map,ierr),*ierr);

  // residual vectors
  PHIST_CHK_IERR(SUBR(mvec_create)(&R,map,numSys,ierr),*ierr);
  // vectors V(:,j) for each of the systems
  PHIST_CHK_IERR(SUBR(mvec_create)(&V,map,numSys,ierr),*ierr);
  // vectors W=A*V(:,j) for each of the systems
  PHIST_CHK_IERR(SUBR(mvec_create)(&W,map,numSys,ierr),*ierr);

  bool anyConverged=false;
  bool anyFailed=false;

  // we exit the while loop as soon as one system converges or reaches its
  // maximum permitted number of iterations. The decision about what to do
  // next is then left to the caller.
  while (1)
    {
    //R=B-A*X;
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one,B,st::zero(),R,ierr), *ierr);
    PHIST_CHK_IERR(Op->apply(-st::one(),Op->A, X, st::one(), R, ierr), *ierr);
  
    PHIST_CHK_IERR(mvec_norm2)(R,beta,ierr),*ierr);
 
    for (int i=0;i<numSys;i++)
      {
      if (S_array[i]->normB_<mt::zero())
        {
        PHIST_OUT(PHIST_ERROR,"state object %d  not initialized correctly",i);
        *ierr=-1;
        return;
        }
      S_array[i]->normR_=beta[i];
      if (S_array[i]->normR_ <= S_array[i]->normB_*tol;
        {
        anyConverged=true;
        S_array[i]->ierr=0;
        }
      if (S_array[i]->numIters_>S_array[i]
        {
        anyFailed=true;
        S_array[i]->ierr=2;
        }
      }
  
    if (anyConverged)
      {
      *ierr=0;
      break;
      }
      
    if (anyFailed)
      {
      *ierr=1;
      break;
      }

  // while loop for inner GMRES iterations until restart
  while (1)
    {
    // Arnoldi - build orthogonal basis V and upper Hessenberg matrix H.
    // jmax is chosen so that the next system reaches maxBas.
    int jmax = S_array[0]->maxBas-S_array[0]->curDimV_;
    for (int i=1;i<numSys;i++)
      {
      jmax = std::min(jmax, S_array[i]->maxBas-S_array[i]->curDimV_);
      }
    for (int j_dum=0; j_dum<jmax; j_dum++)
      {
      //W=A*(M\V(:,j);
      PHIST_CHK_IERR(Op->apply(st::one(),Op->A, V, st::zero(), W, ierr), *ierr);
      
      // Arnoldi update. TODO - we could save some messages by
      // clustering the communication in orthog
      // TODO - maybe we could parallelize this loop by an OpenMP section
      //        so that the small stuff can be done in parallel too?
      for (int i=0;i<numSys;i++)
        {
        TYPE(gmresState) *S = S_array+i;
        int j=S->curDim_;
        //TODO - blocking of vectors in orthog to avoid communication??
        
        // view the existing basis vectors as Vprev
        PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Vprev,0,j-1,ierr),*ierr);
        // view the next V vector as Vj
        PHIST_CHK_IERR(SUBR(mvec_view_block)(S->V_,&Vj,j,j,ierr),*ierr);
        // view H(j,j) as R1
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_,&R1,j,j,j,j,ierr),*ierr);
        // view H(1:j,j) as R2
        PHIST_CHK_IERR(SUBR(sdMat_view_block)(S->H_,&R2,0,j-1,j,j,ierr),*ierr);
        //orthogonalize
        PHIST_CHK_IERR(SUBR(orthog)(Vprev,V,R1,R2,ierr),*ierr);
        
        //    % update QR factorization of H

        // note: it is OK to work on raw data of serial dense matrices because
        // they are typically not modified on an accelerator, however, we must
        // somehow implement a check of this in the kernel interfaces (TODO)
        // note also that we will access H_col[jj], which is strictly speaking
        // R1, but by construction they should be aligned in memory.
        ST *Hcol=NULL;
        lidx_t ldH
        PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R2,&H_col,&ldH,ierr),*ierr);

        //    % apply previous (j-1) transformations to columns j
        for (int jj=0;jj<j-1;jj++)
          {
          ST htmp = S->cs_[jj]*Hcol[jj] + 
                    S->sn_[jj]*Hcol[jj+1]; // H(j-1,j) in the last step, which is R1 above
          Hcol[jj+1] = -st::conj(S->sn_[jj])*Hcol[jj] + S->cs_[jj]*Hcol[jj+1];
          Hcol[jj]   = htmp;
          }

        // new Givens rotation for eliminating H(j+1,j)
        [cs,sn,htmp]=give(Hcol(j),Hcol(j+1));
    trans(1,j)=cs;
    trans(2,j)=sn;
    % eliminate H(j+1,j)
    Hcol(j) = htmp;
    Hcol(j+1)=0;
    % apply to RHS
    tmp=cs*rs(j);
    rs(j+1) = -conj(sn).*rs(j);
    rs(j)=tmp;

    H(1:j+1,j)=Hcol;

    if (debug)
        disp('transformed H');
        H(1:k*(j+1),1:j)
        disp('rs');
      rs(1:j+1,:)
    end
    relres = abs(rs(j+1))./bnorm;
    resvec=[resvec;relres];
    if (verbose)
      print_iter(iter,relres);
    end
    if (max(relres)<=tol || iter>=maxIter)
      %disp('convergence - exit Arnoldi');
      break;
    end
    %disp(sprintf('j=%d, cs=%8.4f sn=%8.4f rs=%8.4f',j,trans(1,j),trans(2,j),rs(j)));
  end % Arnoldi
  y = triu(H(1:j,1:j))\rs(1:j);
  X = X+apply_op(V(:,1:j)*y,M);    
    }//while
  }


