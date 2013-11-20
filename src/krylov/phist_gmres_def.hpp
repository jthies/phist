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
  V->cs_ = new ST[maxBas];
  V->sn_ = new ST[maxBas];
  V->rs_ = new ST[maxBas];
  S->curDim_=0;
  S->curIter_=0;
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
SUBR(gmresState_reset)(TYPE(gmresState_ptr) S, TYPE(const_mvec_ptr) x0,int *ierr)
  {
#include "phist_std_typedefs.hpp"  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  MT rs0;
  TYPE(mvec_ptr) v0 = NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V_,&v0,0,0,ierr),*ierr);
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
        TYPE(gmresState)* array_of_states,
        int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // map defining how to create new vector objects
  const_map_t* map=NULL;
  // multi-vectors to work with several of the systems simultaneously
  TYPE(mvec_ptr) R=NULL, V=NULL, W=NULL;
  // views of the active columns
  TYPE(mvec_ptr) Rv=NULL, Vv=NULL, Wv=NULL,Xv=NULL,Bv=NULL;
  // for the orthog routine
  TYPE(sdMat_ptr) R1=NULL,R2=NULL; 
  
  // number of systems we're iterating on. This value may decrease
  // as systems converge, until at least 'stopIf' of them are
  // converged and we exit the subroutine.
  int numActive;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&numActive,ierr),*ierr);
  
  // figure out how new vectors are built
  PHIST_CHK_IERR(SUBR(mvec_get_map)(X,&map,ierr),*ierr);

  // residual vectors
  PHIST_CHK_IERR(SUBR(mvec_create)(&R,map,numSys,ierr),*ierr);

  bool anyConverged=false;
  bool anyFailed=false;

  // we exit the while loop as soon as one system converges or reaches its
  // maximum permitted number of iterations. The decision about what to do
  // next is then left to the caller.
  while (anyConverged==false && anyFailed==false)
    {
    // view the active columns of R, X and B
    PHIST_CHK_IERR(SUBR(mvec_view_block)(R,&Rv,0,numActive-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(X,&Xv,0,numActive-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(B,&Bv,0,numActive-1,ierr),*ierr);
    
    //R=B-A*X;
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one,Bv,st::zero(),Rv,ierr), *ierr);
    PHIST_CHK_IERR(Op->apply(-st::one(),Op->A, Xv, st::one(), Rv, ierr), *ierr);
  //beta = ||R||_2
  if (max(beta/bnorm)<=tol)
    flag=0;
    break;
  end
  if (iter>maxIter)
    flag=-1;
%    disp('max iter exceeded');
    break;
  end
  // reset gmresState -> normalize initial guess for next maxBas steps
  PHIST_CHK_IERR(SUBR(gmresState_reset)(S,R,ierr),*ierr);
  // Arnoldi - build orthogonal basis V and upper Hessenberg matrix H
  for (int j=1; j<S->maxBas; j++)
    S->curIter_++;
    //W=A*(M\V(:,j);
    // Arnoldi update
//    [Vcol,Hcol]=orthog(V(:,1:j),W,relax);
//    V(:,j+1)=Vcol;
//    % update QR factorization of H

//    % apply previous (j-1) transformations to columns j
    for jj=1:j-1 % apply Givens rotation
      htmp = trans(1,jj)*Hcol(jj) + ...
           trans(2,jj)*Hcol(jj+1);
      Hcol(jj+1) = -conj(trans(2,jj))*Hcol(jj) + trans(1,jj)*Hcol(jj+1);
      Hcol(jj)   = htmp;
    end

    % new Givens rotation for eliminating H(j+1,j)
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
  // reorder the columns of X to reflect the original ordering of B
  }


