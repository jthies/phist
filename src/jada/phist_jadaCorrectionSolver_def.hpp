/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! create a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_create)(TYPE(jadaCorrectionSolver_ptr) *me, phist_jadaOpts opts,
        phist_const_map_ptr map, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
  *me = new TYPE(jadaCorrectionSolver);
  (*me)->method_ = opts.innerSolvType;
  (*me)->leftPrecon=nullptr;
  (*me)->rightPrecon=nullptr;

  int innerSolvBlockSize=opts.innerSolvBlockSize;
  if (innerSolvBlockSize<0) innerSolvBlockSize=opts.blockSize;
  PHIST_CHK_IERR( *iflag = (innerSolvBlockSize <= 0) ? PHIST_INVALID_INPUT : 0, *iflag);

  (*me)->hermitian           = opts.symmetry==phist_HERMITIAN? 1:0;
  #ifndef IS_COMPLEX
  if (opts.symmetry==phist_COMPLEX_SYMMETRIC) (*me)->hermitian = 1;
  #endif
  (*me)->innerSolvBlockSize_ = innerSolvBlockSize;
  (*me)->preconSkewProject=opts.preconSkewProject;
  if (opts.preconOp==nullptr) (*me)->preconSkewProject=0;

  if ((*me)->method_==phist_GMRES||(*me)->method_==phist_MINRES)
  {
    (*me)->blockedGMRESstates_  = new TYPE(blockedGMRESstate_ptr)[(*me)->innerSolvBlockSize_];
    (*me)->innerSolvMaxBas_ = opts.innerSolvMaxBas;
    if ((*me)->innerSolvMaxBas_<0) (*me)->innerSolvMaxBas_=opts.innerSolvMaxIters;
    PHIST_CHK_IERR(SUBR(blockedGMRESstates_create)((*me)->blockedGMRESstates_, innerSolvBlockSize, map, (*me)->innerSolvMaxBas_, iflag), *iflag);
    // note: GMRES and MINRES implement left preconditioning only in phist
    (*me)->leftPrecon=(TYPE(linearOp_ptr))opts.preconOp;
  }
  else if ((*me)->method_==phist_BICGSTAB)
  {
    if (opts.preconOp!=nullptr)
    {
      PHIST_SOUT(PHIST_WARNING, "preconditioning not implemented for QMR and BiCGStab, we suggest to use IDRS instead. Your preconditioner will be ignored.")
    }
  }
  else if ((*me)->method_==phist_QMR)
  {
    (*me)->rightPrecon=(TYPE(linearOp_ptr))opts.preconOp;
  }
  else if ((*me)->method_==phist_IDRS)
  {
    PHIST_SOUT(PHIST_INFO, "TROET: got IDRs as inner solver!");
    (*me)->innerSolvMaxBas_    = opts.innerSolvMaxBas;
    if ((*me)->innerSolvMaxBas_<0) (*me)->innerSolvMaxBas_=4;
    (*me)->rightPrecon=(TYPE(linearOp_ptr))opts.preconOp;
  }
  else if ((*me)->method_==phist_CARP_CG)
  {
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED, *iflag);
  }
  else if ((*me)->method_==phist_USER_LINSOLV)
  {
    if (opts.customSolver_run==nullptr && opts.customSolver_run1==nullptr)
    {
      *iflag=-88;
    }
    typedef void (*customSolver_run_funptr_type)(         void*  customSolverData,
                                    void const*    A_op,       void const*    B_op,
                                    void const*    Qtil,       void const*    BQtil,
                                    const double sigma_r[],    const double sigma_i[],
                                    TYPE(const_mvec_ptr)  res, const int resIndex[],
                                    const double        tol[], int       maxIter,
                                    TYPE(mvec_ptr)        t,
                                    int robust,                int abortAfterFirstConvergedInBlock,
                                    int *iflag);
    typedef void (*customSolver_run1_funptr_type)(        void*  customSolverData,
                                    void const*    A_op,     void const*    B_op,
                                    void const*    Qtil,     void const*    BQtil,
                                    double   sigma_r,        double sigma_i, 
                                     TYPE(const_mvec_ptr)  res,
                                    const double           tol,    int                   maxIter,
                                    TYPE(mvec_ptr)        t,
                                    int robust,
                                    int *                 iflag);
    (*me)->customSolver_=opts.customSolver;
    (*me)->customSolver_run= (customSolver_run_funptr_type) opts.customSolver_run;
    (*me)->customSolver_run1= (customSolver_run1_funptr_type) opts.customSolver_run1;
  }
  else if ((*me)->method_==phist_NO_LINSOLV)
  {
    (*me)->preconSkewProject=opts.preconSkewProject;
    (*me)->rightPrecon=    (*me)->rightPrecon=(TYPE(linearOp_ptr))opts.preconOp;
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "method %d (%s) not implemented\n",(int)(*me)->method_, linSolv2str((*me)->method_));
    *iflag=PHIST_NOT_IMPLEMENTED;
  }
}

//! delete a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_delete)(TYPE(jadaCorrectionSolver_ptr) me, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  if (me->method_==phist_GMRES || me->method_==phist_MINRES)
  {
    PHIST_CHK_IERR(SUBR(blockedGMRESstates_delete)(me->blockedGMRESstates_, me->innerSolvBlockSize_, iflag), *iflag);
    delete[] me->blockedGMRESstates_;
  }
  else if (me->method_==phist_CARP_CG)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
  }
  else if (me->method_==phist_USER_LINSOLV)
  {
  }
  delete me;
}


//! calculate approximate solutions to given set of jacobi-davidson correction equations
//!
//! arguments:
//! jdCorrSolver    the jadaCorrectionSolver object
//! AB_op            matrix A passed to jadaOp_create
//! B_op            matrix B passed to jadaOp_create
//! Qtil            projection vectors V passed to jadaOp_create
//! BQtil           projection vectors BV passed to jadaOp_create
//! sigma           (pos.!) shifts, -sigma[i], i in {1, ..., nvec} is passed to the jadaOp
//! res             JD residuals, e.g. rhs of the correction equations
//! resIndex        if not NULL, specifies permutation of the residual array to avoid unnecessary copying in the jada-algorithm
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! t               returns approximate solution vectors
//! iflag            a value > 0 indicates the number of systems that have not converged to the desired tolerance
void SUBR(jadaCorrectionSolver_run)(TYPE(jadaCorrectionSolver_ptr) me,
                                    TYPE(const_linearOp_ptr)    AB_op,     TYPE(const_linearOp_ptr)    B_op, 
                                    TYPE(const_mvec_ptr)  Qtil,     TYPE(const_mvec_ptr)  BQtil,
                                    const _ST_            sigma[],  TYPE(const_mvec_ptr)  res,      const int resIndex[], 
                                    const _MT_            tol[],    int                   maxIter,
                                    TYPE(mvec_ptr)        t,
                                    int useIMGS,                    int abortAfterFirstConvergedInBlock,
                                    int preconUpdate,               int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  // total number of systems to solve. Note that we may get more vectors in res but only should
  // consider res(resIndex(1:totalNumSys),:)
  int totalNumSys, totalNumRHS, numProj;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(t, &totalNumSys, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(res, &totalNumRHS, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Qtil, &numProj, iflag), *iflag);
  
  // this check may fail: subspacejada expects us to grab the RHS's needed from res(:,resIndex[0:totalNumSys-1])
  //PHIST_CHK_IERR(*iflag=totalNumSys==totalNumRHS?0:PHIST_INVALID_INPUT,*iflag);

  if (me->method_==phist_USER_LINSOLV)
  {
    if (totalNumSys==1 && me->customSolver_run1!=nullptr)
    {
      TYPE(mvec_ptr) res0=(TYPE(mvec_ptr))res;
      if (resIndex)
      {
        res0=nullptr;
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res,&res0,resIndex[0],resIndex[0],iflag),*iflag);
      }
      phist::MvecOwner<_ST_> _res0(res0!=res?res0:nullptr);
      PHIST_CHK_IERR(me->customSolver_run1(me->customSolver_,AB_op,B_op,Qtil,BQtil,(double)st::real(sigma[0]),
        (double)st::imag(sigma[0]), res0,
        (double)tol[0],maxIter,t,useIMGS,iflag),*iflag);
    }
    else if (me->customSolver_run!=nullptr)
    {
      double sr[totalNumSys], si[totalNumSys],dtol[totalNumSys];
      for (int i=0;i<totalNumSys;i++)
      {
        sr[i]=st::real(sigma[i]);
        si[i]=st::imag(sigma[i]);
        dtol[i]=(double)tol[i];
      }
      PHIST_CHK_IERR(me->customSolver_run(me->customSolver_,AB_op,B_op,Qtil,BQtil,sr,si,res,resIndex,
        dtol,maxIter,t,useIMGS,abortAfterFirstConvergedInBlock,iflag),*iflag);
    }
    else if (me->customSolver_run1!=nullptr)
    {
      static bool first_time=true;
      if (first_time)
      {
        PHIST_SOUT(PHIST_WARNING,"custom solver requested for block size %d. As you provided only\n"
                                 "a function for a single system, I will use it one-by-one.\n",totalNumSys);
        first_time=false;
      }
      TYPE(mvec_ptr) resj=nullptr,tj=nullptr;
      for (int j=0; j<totalNumSys; j++)
      {
        int jres = (resIndex==0)? j: resIndex[j];
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res,&resj,jres,jres,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))t,&tj,j,j,iflag),*iflag);
        
        PHIST_CHK_IERR(me->customSolver_run1(me->customSolver_,AB_op,B_op,Qtil,BQtil,(double)st::real(sigma[j]),
        (double)st::imag(sigma[j]), resj, (double)tol[0],maxIter,tj,useIMGS,iflag),*iflag);
      }
      if (resj!=nullptr) PHIST_CHK_IERR(SUBR(mvec_delete)(resj,iflag),*iflag);
      if (  tj!=nullptr) PHIST_CHK_IERR(SUBR(mvec_delete)(  tj,iflag),*iflag);
    }
    else
    {
      PHIST_SOUT(PHIST_ERROR,"custom solver requested but function not set in jadaOpts struct\n");
      *iflag=PHIST_BAD_CAST;
    }
    return;
  }

  // shifts to pass to the preconditioner, for the moment just set them to 0,
  // if we want to use preconditioners that can handle changing shifts we have
  // to set these depending on resIndex etc. below
  std::vector<_ST_> preconShifts(std::max(totalNumSys,numProj), st::zero());

  // if there is a preconditioner and the options indicate that we want to apply it
  // as a projected operator, create a copy of the last block of Q/BQ and compute  
  // P\Q to be used for the projection. In theory it may be better to project out  
  // all of the vectors in Q also from the preconditioner, but in practice this means
  // we have to store P\Q.
  TYPE(mvec_ptr) q=nullptr, Bq=nullptr;
  // if we do a manual projection after the preconditioner (skew-projection) or
  // are asked to update the preconditioner we set extra pointers q,Bq. In the 
  // case of preconditioner updates it may or may not be that the preconditioner
  // uses the projection space q,Bq.
  if (me->preconSkewProject!=0 || preconUpdate)
  {
    int imax=numProj-1;
    int imin=std::max(0,imax-totalNumSys);
    PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))Qtil,&q,imin,imax,iflag),*iflag);
    if (BQtil!=nullptr) {PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))BQtil,&Bq,imin,imax,iflag),*iflag);}
    else             Bq=q;
  }

  // make sure these vectors get deleted at the end of the scope
  phist::MvecOwner<_ST_> _q(q!=Qtil?q:nullptr),_Bq(( (Bq!=q)&&(Bq!=BQtil) )?Bq:nullptr);
  
  if (preconUpdate!=0)
  {
    // select a single shift for the preconditioner (various strategies could be used)
    _ST_ prec_sigma=sigma[0];
    if (me->leftPrecon!=nullptr)
    {
      // not all preconditioners may support this, so only warn on non-zero output
      SUBR(precon_update)(me->leftPrecon,prec_sigma,q,Bq,iflag);
      if (*iflag) PHIST_SOUT(PHIST_WARNING,"precon_update returned non-zero code %d\n"
                                           "(file %s, line %d). Ignoring return value.\n",
                                           *iflag,__FILE__,__LINE__);
    }
  }

  // if no Krylov method is used, only apply the preconditioner. This is also known as Olsen's method.
  if (me->method_==phist_NO_LINSOLV ||
      me->method_==phist_QMR        ||
      me->method_==phist_BICGSTAB   ||
      me->method_==phist_IDRS)
  {
    TYPE(mvec_ptr) _t=t,_res=(TYPE(mvec_ptr))res;
    int jlower=0,jupper=totalNumSys-1;
    if (resIndex!=nullptr)
    {
      // only implemented up to now for contiguous and increasing resIndex
      bool valid_resIndex=true;
      for (int i=1; i<totalNumSys; i++) valid_resIndex|=(resIndex[i]==resIndex[i-1]+1);
      PHIST_CHK_IERR(*iflag=valid_resIndex?0:PHIST_NOT_IMPLEMENTED,*iflag);
      jlower=resIndex[jlower];
      jupper=resIndex[jupper];
      _res=nullptr;
      PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res,&_res,jlower,jupper,iflag),*iflag);
    }
    // make sure the view gets deleted at the end of the scope
    phist::MvecOwner<_ST_> __res(_res!=res?_res:nullptr);
    // wrap the preconditioner so that apply_shifted is called
    TYPE(linearOp) jadaPrecR, *jadaPrecRight=nullptr;

    if (me->leftPrecon!=nullptr)
    {
      PHIST_SOUT(PHIST_ERROR,"left preconditioning only implemented for GMRES and MINRES in jadaCorrectionSolver right now.\n");
      PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
    }

    if (me->rightPrecon!=nullptr)
    {
      // create preconditioner with totalNumSys zero shifts because we need to apply it to the RHS as well and
      // we check internally that sufficient shifts are given.

      // shifts to pass to the preconditioner, for the moment just set them to 0,
      // if we want to use preconditioners that can handle changing shifts we have
      // to set these depending on resIndex etc. below
      std::vector<_ST_> preconShifts(std::max(totalNumSys,numProj), st::zero());

      PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->rightPrecon,q,Bq,&preconShifts[0],totalNumRHS,&jadaPrecR,me->preconSkewProject,iflag),*iflag);
    }


      if (me->method_==phist_NO_LINSOLV)
      {
        if (jadaPrecRight!=nullptr)
        {
          // apply this preconditioner
          PHIST_CHK_IERR(jadaPrecRight->apply(st::one(),jadaPrecRight->A, _res, st::zero(), _t,iflag),*iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),res,st::zero(),t,iflag),*iflag);
        }

      }
      else if (me->method_==phist_QMR || me->method_==phist_BICGSTAB || me->method_==phist_IDRS)
      {
        // make sure the inner and outer block sizes are the same, and there is no permutation of the residual vectors.
        // For GMRES and MINRES below we don't have this restriction
        int k = std::min(me->innerSolvBlockSize_, totalNumSys);
        PHIST_CHK_IERR(*iflag=(k==totalNumSys)?0:PHIST_NOT_IMPLEMENTED,*iflag);
        // note: we already checked the permutation of the res (rhs vectors) above is consecutive and
        // extracted _res as a view of those columns. Check that _t has the same number of columns.
        int nct;
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(_t,&nct,iflag),*iflag);
        PHIST_CHK_IERR(*iflag=(nct==totalNumSys)?0:PHIST_INVALID_INPUT,*iflag);

        // we need a jadaOp
        TYPE(linearOp) jadaOp;
        PHIST_CHK_IERR(SUBR(jadaOp_create)(AB_op, B_op, Qtil, BQtil, &sigma[0], k, &jadaOp, iflag), *iflag);
        int nIter=maxIter;
        // pass in the last k columns of Qtil (the current approximate eigenspace we're solving for)
        TYPE(mvec_ptr) V=nullptr;
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))Qtil,&V,std::max(numProj-k,0),numProj-1,iflag),*iflag);
        phist::MvecOwner<_ST_> _V(V);
        if (me->method_==phist_IDRS)
        {
          int s = me->innerSolvMaxBas_;
          PHIST_SOUT(PHIST_INFO, "run IDR(%d) for maxIter=%d iterations and %d RHS\n", s, nIter, k);//TROET
          MT min_tol=tol[0]; for (int i=1; i<k; i++) min_tol=std::min(min_tol,tol[i]);
          PHIST_CHK_NEG_IERR(SUBR(blockedIDRs_iterate)(&jadaOp, jadaPrecRight, _res,_t, V, k, &nIter, min_tol, s, iflag),*iflag);
        }
        else if (me->method_==phist_BICGSTAB)
        {
          PHIST_CHK_NEG_IERR(SUBR(blockedBiCGStab_iterate)(&jadaOp, jadaPrecRight, _res,_t, V, k, &nIter, tol, iflag),*iflag);
        }
        else if (me->method_==phist_QMR)
        {
          int sym=me->hermitian;
          PHIST_CHK_NEG_IERR(SUBR(blockedQMR_iterate)(&jadaOp, jadaPrecRight, _res,_t, V, k, &nIter, tol, sym, iflag),*iflag);
        }
      }

      if (jadaPrecRight!=nullptr)
      {
        // delete the preconditioner (wrapper) again
        PHIST_CHK_IERR(jadaPrecRight->destroy(jadaPrecRight,iflag),*iflag);
      }
      *iflag=0;
      return;
    }

    /*
    else: GMRES or MINRES, our default approach.
    This implementation is more involved, it allows solving
    a number of linear systems by dynamically swapping converged
    ones for unconverged ones. This could be useful if a rather
    large block size is used for the outer JD solver, but in practice
    this hasn't been done a lot.
    */

  // MINRES or GMRES implemented below, these are more fancy as they allow restarting,
  // swapping out converged systems etc. In the future we may abolish all of that and treat
  // them in the same way as QMR above. CARP-CG is not implemented at all right now.

  PHIST_CHK_IERR(*iflag = (maxIter <= 0) ? -1 : 0, *iflag);

  // set solution vectors to zero to add them up later
  PHIST_CHK_IERR(SUBR(mvec_put_value)(t, st::zero(), iflag), *iflag);

  // make sure all states are reset
  for(int i = 0; i < me->innerSolvBlockSize_; i++)
  {
    PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(me->blockedGMRESstates_[i], nullptr, nullptr, iflag), *iflag);
  }

  // current and maximal block dimension
  int max_k = me->innerSolvBlockSize_;
  int k = max_k;

  // index of currently iterated systems in all systems to solve
  std::vector<int> index(max_k);

  // shifts currently in use
  std::vector<_ST_> currShifts(max_k, st::zero());

  // we need a jadaOp
  TYPE(linearOp) jadaOp;
  PHIST_CHK_IERR(SUBR(jadaOp_create)(AB_op, B_op, Qtil, BQtil, &currShifts[0], k, &jadaOp, iflag), *iflag);

    TYPE(linearOp) jadaPrecL, *jadaPrecLeft=nullptr;

  // the RHS used, if there is no left preconditioning this is just a pointer to res, otherwise
  // it is allocated and the preconditioner is applied
  TYPE(mvec_ptr) rhs = (TYPE(mvec_ptr))res;

  // wrap the preconditioner so that apply_shifted is called
  if (me->leftPrecon!=nullptr)
  {
    // create preconditioner with totalNumSys zero shifts because we need to apply it to the RHS as well and
    // we check internally that sufficient shifts are given.
    PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->leftPrecon,q,Bq,&preconShifts[0],totalNumRHS,&jadaPrecL,me->preconSkewProject,iflag),*iflag);
    rhs=nullptr;
    PHIST_CHK_IERR(SUBR(mvec_clone_shape)(&rhs,res,iflag),*iflag);
    PHIST_CHK_IERR(jadaPrecL.apply(st::one(),jadaPrecL.A,res,st::zero(),rhs,iflag),*iflag);
    jadaPrecLeft=&jadaPrecL;
    PHIST_CHK_IERR(SUBR(jadaOp_set_leftPrecond)(&jadaOp,jadaPrecLeft,iflag),*iflag);
  }
  
  // in case a preconditioned copy was made, delete it at the end of the scope
  phist::MvecOwner<_ST_> _rhs(rhs!=res?rhs:nullptr);

  int nextSystem = 0;

  // just used to prevent inifite loops
  int finiteLoopCounter = 0;

  // we need some views
  TYPE(mvec_ptr) res_i = nullptr;
  TYPE(mvec_ptr) t_i   = nullptr;

  // number of unconverged systems to be returned
  int nUnconvergedSystems = 0;

// put all iterations in one big compute task; this speeds up the tests with ghost (significantly)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  // iterate while there's at least one system remaining
  while( k > 0 )
  {
    // prevent bugs from causing inifite looping
    if( finiteLoopCounter++ > maxIter )
      break;

    // find unoccupied states
    for(int i = 0; i < max_k; i++)
    {
      if( nextSystem >= totalNumSys )
        break;
      if( me->blockedGMRESstates_[i]->status == -2 )
      {
        // setup the next system waiting to be solved
        int ind = nextSystem;
        if( resIndex != nullptr ) ind = resIndex[ind];
        PHIST_CHK_IERR(SUBR(mvec_view_block)(rhs, &res_i, ind, ind, iflag), *iflag);
        PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(me->blockedGMRESstates_[i], res_i, nullptr, iflag), *iflag);

        me->blockedGMRESstates_[i]->tol = tol[nextSystem];
        index[me->blockedGMRESstates_[i]->id] = nextSystem;

        nextSystem++;
      }
    }

    // gather systems that are waiting to be iterated
    std::vector<TYPE(blockedGMRESstate_ptr)> activeStates;
    int firstId = max_k;
    if (max_k>0) PHIST_SOUT(PHIST_VERBOSE, "Iterating systems:");
    for(int i = 0; i < max_k; i++)
    {
      if( std::abs(me->blockedGMRESstates_[i]->status) == 1 )
      {
#ifdef IS_COMPLEX
        PHIST_SOUT(PHIST_VERBOSE, "\t%d (%f%+fi)", index[me->blockedGMRESstates_[i]->id], 
        -st::real(sigma[index[me->blockedGMRESstates_[i]->id]]),
        -st::imag(sigma[index[me->blockedGMRESstates_[i]->id]]));
#else
        PHIST_SOUT(PHIST_VERBOSE, "\t%d (%f)", index[me->blockedGMRESstates_[i]->id], -sigma[index[me->blockedGMRESstates_[i]->id]]);
#endif
        activeStates.push_back(me->blockedGMRESstates_[i]);
        firstId = std::min(firstId,activeStates.back()->id);
      }
    }
    if (max_k>0) PHIST_SOUT(PHIST_VERBOSE, "\n");
    k = activeStates.size();

    // set correct shifts
    for(int i = 0; i < k; i++)
    {
      currShifts[activeStates[i]->id - firstId] = -sigma[index[activeStates[i]->id]];
    }

    // actually iterate
    if( me->method_==phist_MINRES )
    {
      int nIter=maxIter;
      PHIST_CHK_NEG_IERR(SUBR(blockedMINRESstates_iterate)(&jadaOp, nullptr, &activeStates[0], k, &nIter, iflag), 
      *iflag);
    }
    else if (me->method_==phist_CARP_CG)
    {
      *iflag=PHIST_NOT_IMPLEMENTED;
    }
    else
    {
      int nIter=maxIter;
      PHIST_CHK_NEG_IERR(SUBR(blockedGMRESstates_iterate)(&jadaOp, nullptr,&activeStates[0], k, &nIter, useIMGS, iflag), *iflag);
    }

    // optimization to use always full blocks and ignore tolerances
    if( abortAfterFirstConvergedInBlock && k > 0 )
    {
      int ind = index[activeStates[0]->id];
      PHIST_CHK_IERR(SUBR(mvec_view_block)(t, &t_i, ind, ind+k-1, iflag), *iflag);
      _MT_ tmp[k];
      PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(&activeStates[0], k, nullptr, t_i, tmp, false, iflag), 
      *iflag);
      for(int i = 0; i < k; i++)
      {
        PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(activeStates[i], nullptr, nullptr, iflag), *iflag);
      }
    }
    else
    {
      // check the status of the systems
      for(int i = 0; i < k; i++)
      {
        if( activeStates[i]->status == 0 || activeStates[i]->status >= 2 )
        {
          // update solution
          int ind = index[activeStates[i]->id];
          PHIST_CHK_IERR(SUBR(mvec_view_block)(t, &t_i, ind, ind, iflag), *iflag);
          _MT_ tmp;
          PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(&activeStates[i],1,nullptr, t_i, &tmp, false, iflag), 
          *iflag);

          if( activeStates[i]->status == 3 )
          {
            nUnconvergedSystems++;
          }
          else if( activeStates[i]->status == 2 )
          {
            // prepare restart
            PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(activeStates[i], nullptr, t_i, iflag), *iflag);
            continue;
          }
/*
#ifdef PHIST_TESTING
{
  // determine real residual for comparison
  currShifts[0] = -sigma[ind];
  TYPE(mvec_ptr) res_i = nullptr;
  int resInd = ind;
  if( resIndex != nullptr )
    resInd = resIndex[resInd];
  PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res, &res_i, resInd, resInd, iflag), *iflag);
  _MT_ nrm0;
  PHIST_CHK_IERR(SUBR(mvec_norm2)(res_i, &nrm0, iflag), *iflag);
  PHIST_CHK_IERR(jadaOp.apply(-st::one(), jadaOp.A, t_i, st::one(), res_i, iflag), *iflag);
  _MT_ nrm;
  PHIST_CHK_IERR(SUBR(mvec_norm2)(res_i, &nrm, iflag), *iflag);
  PHIST_SOUT(PHIST_INFO,"est. / exp. residual of system %d: %8.4e / %8.4e\n", ind, tmp, nrm/nrm0);
}
#endif
*/

          // reset to be free in the next iteration
          PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(activeStates[i], nullptr, nullptr, iflag), *iflag);
        }
      }
    }

  }


  // normalize result vectors, TODO: should be done in updateSol/pgmres?
  _MT_ tmp[totalNumSys];
  PHIST_CHK_IERR(SUBR(mvec_normalize)(t, tmp, iflag), *iflag);

  // delete views
  PHIST_CHK_IERR(SUBR(mvec_delete)(res_i, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(t_i,   iflag), *iflag);
PHIST_TASK_END(iflag)

  // delete the jadaOp and preconditioner wrappers
  if (jadaPrecLeft!=nullptr)
  {
    PHIST_CHK_IERR(jadaPrecLeft->destroy(jadaPrecLeft, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR(jadaOp_delete)(&jadaOp, iflag), *iflag);
  *iflag = nUnconvergedSystems;
}
