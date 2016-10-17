//! create a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_create)(TYPE(jadaCorrectionSolver_ptr) *me, phist_jadaOpts opts,
        phist_const_map_ptr map, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
  *me = new TYPE(jadaCorrectionSolver);
  (*me)->method_ = opts.innerSolvType;
  if ((*me)->method_==phist_GMRES||(*me)->method_==phist_MINRES)
  {
    PHIST_CHK_IERR( *iflag = (opts.innerSolvBlockSize <= 0) ? -1 : 0, *iflag);

    (*me)->gmresBlockDim_ = opts.innerSolvBlockSize;
    (*me)->blockedGMRESstates_  = new TYPE(blockedGMRESstate_ptr)[(*me)->gmresBlockDim_];
    PHIST_CHK_IERR(SUBR(blockedGMRESstates_create)((*me)->blockedGMRESstates_, opts.innerSolvBlockSize, map, opts.innerSolvMaxBas, iflag), *iflag);
    (*me)->rightPrecon=(TYPE(const_linearOp_ptr))opts.preconOp;
    (*me)->preconSkewProject=opts.preconSkewProject;
  }
  else if ((*me)->method_==phist_CARP_CG)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
  }
  else if ((*me)->method_==phist_USER_DEFINED)
  {
    if (opts.customSolver_run==NULL && opts.customSolver_run1==NULL)
    {
      *iflag=-88;
    }
    (*me)->customSolver_=opts.customSolver;
    (*me)->customSolver_run=opts.customSolver_run;
    (*me)->customSolver_run1=opts.customSolver_run1;
  }
  else if ((*me)->method_==phist_NO_LINSOLV)
  {
    (*me)->rightPrecon=(TYPE(const_linearOp_ptr))opts.preconOp;
    (*me)->preconSkewProject=opts.preconSkewProject;
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
    PHIST_CHK_IERR(SUBR(blockedGMRESstates_delete)(me->blockedGMRESstates_, me->gmresBlockDim_, iflag), *iflag);
    delete[] me->blockedGMRESstates_;
  }
  else if (me->method_==phist_CARP_CG)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
  }
  else if (me->method_==phist_USER_DEFINED)
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
                                    bool useIMGS,                   bool abortAfterFirstConvergedInBlock,
                                    int *                 iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  // total number of systems to solve
  int totalNumSys;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(t, &totalNumSys, iflag), *iflag);

  // if there is a preconditioner and the options indicate that we want to apply it
  // as a projected operator, create a copy of the last block of Q/BQ and compute  
  // P\Q to be used for the projection. In theory it may be better to project out  
  // all of the vectors in Q also from the preconditioner, but in practice this means
  // we have to store P\Q.
  TYPE(mvec_ptr) q=NULL, Bq=NULL;
  if (me->rightPrecon && me->preconSkewProject)
  {
    phist_const_map_ptr map;
    int nq;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Qtil,&nq,iflag),*iflag);
    int nqp=totalNumSys;
    int nq0=std::max(0,nq-nqp);
    PHIST_CHK_IERR(SUBR(mvec_get_map)(Qtil,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&q,map,nqp,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_get_block)(Qtil,q,nq0,nq-1,iflag),*iflag);
    if (B_op!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_create)(&Bq,map,nqp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_get_block)(BQtil,Bq,nq0,nq-1,iflag),*iflag);      
    }
  }
  // make sure these vectors get deleted at the end of the scope
  MvecOwner<_ST_> _q(q),_Bq(Bq);
  
  if (me->method_==phist_NO_LINSOLV)
  {
    // wrap the preconditioner so that apply_shifted is called
    if (me->rightPrecon!=NULL)
    {
      std::vector<_ST_> shifts(totalNumSys, st::zero());
      for (int i=0; i<totalNumSys; i++) shifts[i]=-sigma[i];
      if (me->preconSkewProject==0)
      {
        PHIST_CHK_IERR(me->rightPrecon->apply_shifted(st::one(),me->rightPrecon->A, 
                &shifts[0],res,st::one(),t,iflag),*iflag);
      }
      else
      {
        static int first_time=1;
        if (first_time)
        {
          PHIST_SOUT(PHIST_WARNING,"NO_LINSOLV/Olsens method (apply projected preconditioner only)\n"
                                   "requires two preconditioner applications per correction vector\n"
                                   "and may therefore be quite expensive...\n");
          first_time=0;
        }
        TYPE(linearOp) jadaPrec;
        PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->rightPrecon,q,Bq,&shifts[0],totalNumSys,&jadaPrec,iflag),*iflag);
        // apply this preconditioner
        PHIST_CHK_IERR(jadaPrec.apply(st::one(),jadaPrec.A, res, st::zero(), t,iflag),*iflag);
        // delete the preconditioner (wrapper) again
        PHIST_CHK_IERR(jadaPrec.destroy(&jadaPrec,iflag),*iflag);
      }
    }
    else
    {
      static int first_time=1;
      if (first_time)
      {
        PHIST_SOUT(PHIST_WARNING,"NO_LINSOLV/Olsens method (apply projected preconditioner only)\n"
                                 "was chosen in the jadaOpts, but no preconditioner is given. This\n"
                                 "means that the correction equation is not solved!\n");
        first_time=0;
      }
    }
    *iflag=0;
    return;
  }
  if (me->method_==phist_USER_DEFINED)
  {
    int numSys;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(t,&numSys,iflag),*iflag);
    if (numSys==1 && me->customSolver_run1!=NULL)
    {
      PHIST_CHK_IERR(me->customSolver_run1(me->customSolver_,AB_op,B_op,Qtil,BQtil,(double)st::real(sigma[0]),
        (double)st::imag(sigma[0]), res,
        (double)tol[0],maxIter,t,useIMGS,iflag),*iflag);
    }
    else if (me->customSolver_run!=NULL)
    {
      double sr[numSys], si[numSys],dtol[numSys];
      for (int i=0;i<numSys;i++)
      {
        sr[i]=st::real(sigma[i]);
        si[i]=st::imag(sigma[i]);
        dtol[i]=(double)tol[i];
      }
      PHIST_CHK_IERR(me->customSolver_run(me->customSolver_,AB_op,B_op,Qtil,BQtil,sr,si,res,resIndex,
        dtol,maxIter,t,useIMGS,abortAfterFirstConvergedInBlock,iflag),*iflag);
    }
    else
    {
      PHIST_SOUT(PHIST_ERROR,"custom solver requested but function not set in jadaOpts struct\n");
      *iflag=-88;
    }
    return;
  }

  // all other cases implemented below: MINRES, GMRES, CARP-CG

  PHIST_CHK_IERR(*iflag = (maxIter <= 0) ? -1 : 0, *iflag);

  // set solution vectors to zero to add them up later
  PHIST_CHK_IERR(SUBR(mvec_put_value)(t, st::zero(), iflag), *iflag);

  // make sure all states are reset
  for(int i = 0; i < me->gmresBlockDim_; i++)
  {
    PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(me->blockedGMRESstates_[i], NULL, NULL, iflag), *iflag);
  }

  // current and maximal block dimension
  int max_k = me->gmresBlockDim_;
  int k = max_k;

  // index of currently iterated systems in all systems to solve
  std::vector<int> index(max_k);

  // shifts currently in use
  std::vector<_ST_> currShifts(max_k, st::zero());

  // we need a jadaOp
  TYPE(linearOp) jadaOp;
  PHIST_CHK_IERR(SUBR(jadaOp_create)(AB_op, B_op, Qtil, BQtil, &currShifts[0], k, &jadaOp, iflag), *iflag);

    TYPE(linearOp) jadaPrec, *jadaPrecPtr=NULL;

  // wrap the preconditioner so that apply_shifted is called
  if (me->rightPrecon!=NULL)
  {
    if (me->preconSkewProject==0)
    {
      PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->rightPrecon,NULL,NULL,&currShifts[0],k,&jadaPrec,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->rightPrecon,q,Bq,&currShifts[0],k,&jadaPrec,iflag),*iflag);
    }
    jadaPrecPtr=&jadaPrec;
  }

  int nextSystem = 0;

  // just used to prevent inifite loops
  int finiteLoopCounter = 0;

  // we need some views
  TYPE(mvec_ptr) res_i = NULL;
  TYPE(mvec_ptr) t_i   = NULL;

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
        if( resIndex != NULL )
          ind = resIndex[ind];
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res, &res_i, ind, ind, iflag), *iflag);
        PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(me->blockedGMRESstates_[i], res_i, NULL, iflag), *iflag);
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
      PHIST_CHK_NEG_IERR(SUBR(blockedMINRESstates_iterate)(&jadaOp, jadaPrecPtr, &activeStates[0], k, &nIter, iflag), *iflag);
    }
    else if (me->method_==phist_CARP_CG)
    {
      *iflag=PHIST_NOT_IMPLEMENTED;
    }
    else
    {
      int nIter=maxIter;
      PHIST_CHK_NEG_IERR(SUBR(blockedGMRESstates_iterate)(&jadaOp, jadaPrecPtr,&activeStates[0], k, &nIter, useIMGS, iflag), *iflag);
    }

    // optimization to use always full blocks and ignore tolerances
    if( abortAfterFirstConvergedInBlock && k > 0 )
    {
      int ind = index[activeStates[0]->id];
      PHIST_CHK_IERR(SUBR(mvec_view_block)(t, &t_i, ind, ind+k-1, iflag), *iflag);
      _MT_ tmp[k];
      PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(&activeStates[0], k, jadaPrecPtr, t_i, tmp, false, iflag), *iflag);
      for(int i = 0; i < k; i++)
      {
        PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(activeStates[i], NULL, NULL, iflag), *iflag);
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
          PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(&activeStates[i],1,jadaPrecPtr, t_i, &tmp, false, iflag), *iflag);

          if( activeStates[i]->status == 3 )
          {
            nUnconvergedSystems++;
          }
          else if( activeStates[i]->status == 2 )
          {
            // prepare restart
            PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(activeStates[i], NULL, t_i, iflag), *iflag);
            continue;
          }
/*
#ifdef PHIST_TESTING
{
  // determine real residual for comparison
  currShifts[0] = -sigma[ind];
  TYPE(mvec_ptr) res_i = NULL;
  int resInd = ind;
  if( resIndex != NULL )
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
          PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(activeStates[i], NULL, NULL, iflag), *iflag);
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
  // delete the jadaOp and preconditioner wrapper
  PHIST_CHK_IERR(SUBR(jadaOp_delete)(&jadaOp, iflag), *iflag);
  if (jadaPrecPtr!=NULL)
  {
    PHIST_CHK_IERR(jadaPrec.destroy(&jadaPrec, iflag), *iflag);
  }
  *iflag = nUnconvergedSystems;
}

