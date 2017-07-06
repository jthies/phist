/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
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
  (*me)->leftPrecon=NULL; 
  (*me)->rightPrecon=NULL; 

  int innerSolvBlockSize=opts.innerSolvBlockSize;
  if (innerSolvBlockSize<0) innerSolvBlockSize=opts.blockSize;
  PHIST_CHK_IERR( *iflag = (innerSolvBlockSize <= 0) ? PHIST_INVALID_INPUT : 0, *iflag);

  (*me)->innerSolvBlockSize_ = innerSolvBlockSize;

  if ((*me)->method_==phist_GMRES||(*me)->method_==phist_MINRES)
  {
    (*me)->blockedGMRESstates_  = new TYPE(blockedGMRESstate_ptr)[(*me)->innerSolvBlockSize_];
    int innerSolvMaxBas = opts.innerSolvMaxBas;
    if (innerSolvMaxBas<0) innerSolvMaxBas=opts.innerSolvMaxIters;
    PHIST_CHK_IERR(SUBR(blockedGMRESstates_create)((*me)->blockedGMRESstates_, innerSolvBlockSize, map, innerSolvMaxBas, iflag), *iflag);
    (*me)->leftPrecon=(TYPE(linearOp_ptr))opts.preconOp;
    (*me)->preconSkewProject=opts.preconSkewProject;
  }
  else if ((*me)->method_==phist_QMR)
  {
    //TODO: decide wether to use left or right preconditioning with QMR
    (*me)->rightPrecon=(TYPE(linearOp_ptr))opts.preconOp;
    (*me)->preconSkewProject=opts.preconSkewProject;
  }
  else if ((*me)->method_==phist_CARP_CG)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }
  else if ((*me)->method_==phist_USER_LINSOLV)
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
    (*me)->rightPrecon=(TYPE(linearOp_ptr))opts.preconOp;
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
    if (totalNumSys==1 && me->customSolver_run1!=NULL)
    {
      TYPE(mvec_ptr) res0=(TYPE(mvec_ptr))res;
      if (resIndex)
      {
        res0=NULL;
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res,&res0,resIndex[0],resIndex[0],iflag),*iflag);
      }
      MvecOwner<_ST_> _res0(res0!=res?res0:NULL);
      PHIST_CHK_IERR(me->customSolver_run1(me->customSolver_,AB_op,B_op,Qtil,BQtil,(double)st::real(sigma[0]),
        (double)st::imag(sigma[0]), res0,
        (double)tol[0],maxIter,t,useIMGS,iflag),*iflag);
    }
    else if (me->customSolver_run!=NULL)
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
  TYPE(mvec_ptr) q=NULL, Bq=NULL;
  // if we do a manual projection after the preconditioner (skew-projection) or
  // are asked to update the preconditioner we set extra pointers q,Bq. In the 
  // case of preconditioner updates it may or may not be that the preconditioner
  // uses the projection space q,Bq.
  if (me->preconSkewProject!=0 || preconUpdate)
  {
    q=(TYPE(mvec_ptr))Qtil; Bq=(TYPE(mvec_ptr))BQtil;
  }

  if (me->rightPrecon && me->preconSkewProject)
  {
    q=NULL; Bq=NULL;
    phist_const_map_ptr map;
    //int nqp=totalNumSys;
    int nqp=numProj;
    int nq=numProj;
    int nq0=std::max(0,nq-nqp);
    PHIST_CHK_IERR(SUBR(mvec_get_map)(Qtil,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&q,map,nqp,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_get_block)(Qtil,q,nq0,nq-1,iflag),*iflag);
    if (B_op!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_create)(&Bq,map,nqp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_get_block)(BQtil,Bq,nq0,nq-1,iflag),*iflag);      
    }
    else
    {
      Bq=q;
    }
  }

  // make sure these vectors get deleted at the end of the scope
  MvecOwner<_ST_> _q(q!=Qtil?q:NULL),_Bq(( (Bq!=q)&&(Bq!=BQtil) )?Bq:NULL);
  
  if (preconUpdate!=0)
  {
    // select a single shift for the preconditioner (various strategies could be used)
    _ST_ prec_sigma=sigma[0];
    if (me->leftPrecon!=NULL)
    {
      // not all preconditioners may support this, so only warn on non-zero output
      SUBR(precon_update)(me->leftPrecon,prec_sigma,q,Bq,iflag);
      if (*iflag) PHIST_SOUT(PHIST_WARNING,"precon_update returned non-zero code %d\n"
                                           "(file %s, line %d). Ignoring return value.\n",
                                           *iflag,__FILE__,__LINE__);
    }
    if (me->rightPrecon!=NULL)
    {
      // not all preconditioners may support this, so only warn on non-zero output
      SUBR(precon_update)(me->rightPrecon,prec_sigma,q,Bq,iflag);
      if (*iflag) PHIST_SOUT(PHIST_WARNING,"precon_update returned non-zero code %d\n"
                                           "(file %s, line %d). Ignoring return value.\n",
                                           *iflag,__FILE__,__LINE__);
    }
  }
 
  // if no Krylov method is used, only apply the preconditioner. This is also known as Olsen's method.
  if (me->method_==phist_NO_LINSOLV || me->method_==phist_QMR)
  {
    TYPE(mvec_ptr) _t=t,_res=(TYPE(mvec_ptr))res;
    int jlower=0,jupper=totalNumSys-1;
    if (resIndex!=NULL)
    {
      // only implemented up to now for contiguous and increasing resIndex
      bool valid_resIndex=true;
      for (int i=1; i<totalNumSys; i++) valid_resIndex|=(resIndex[i]==resIndex[i-1]+1);
      PHIST_CHK_IERR(*iflag=valid_resIndex?0:PHIST_NOT_IMPLEMENTED,*iflag);
      jlower=resIndex[jlower];
      jupper=resIndex[jupper];
      _res=NULL;
      PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res,&_res,jlower,jupper,iflag),*iflag);
    }
    // make sure the view gets deleted at the end of the scope
    MvecOwner<_ST_> __res(_res!=res?_res:NULL);
    // wrap the preconditioner so that apply_shifted is called
      TYPE(linearOp_ptr) jadaPrec=NULL;

    if (me->rightPrecon!=NULL)
    {
      std::vector<_ST_> shifts(totalNumSys, st::zero());
      if (me->method_==phist_NO_LINSOLV)
      {
        // use actual shifts instead of 0's
        for (int i=0; i<totalNumSys; i++) shifts[i]=-sigma[i];
      }

      PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->rightPrecon,q,Bq,&shifts[0],totalNumSys,jadaPrec,me->preconSkewProject,iflag),*iflag);

    }

      if (me->method_==phist_NO_LINSOLV)
      {
        if (jadaPrec!=NULL)
        {
          // apply this preconditioner
          PHIST_CHK_IERR(jadaPrec->apply(st::one(),jadaPrec->A, _res, st::zero(), _t,iflag),*iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),res,st::zero(),t,iflag),*iflag);
        }

      }
      else if (me->method_==phist_QMR)
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
        int sym=0;
        PHIST_CHK_NEG_IERR(SUBR(blockedQMR_iterate)(&jadaOp, jadaPrec, _res,_t, k, &nIter, tol, sym, iflag),*iflag);
      }

      if (jadaPrec!=NULL)
      {
        // delete the preconditioner (wrapper) again
        PHIST_CHK_IERR(jadaPrec->destroy(jadaPrec,iflag),*iflag);
      }
      *iflag=0;
      return;
    }

  
  // MINRES or GMRES implemented below, these are more fancy as they allow restarting,
  // swapping out converged systems etc. In the future we may abolish all of that and treat
  // them in the same way as QMR above. CARP-CG is not implemented at all right now.

  PHIST_CHK_IERR(*iflag = (maxIter <= 0) ? -1 : 0, *iflag);

  // set solution vectors to zero to add them up later
  PHIST_CHK_IERR(SUBR(mvec_put_value)(t, st::zero(), iflag), *iflag);

  // make sure all states are reset
  for(int i = 0; i < me->innerSolvBlockSize_; i++)
  {
    PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(me->blockedGMRESstates_[i], NULL, NULL, iflag), *iflag);
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

    TYPE(linearOp) jadaPrecL, jadaPrecR, *jadaPrecLeft=NULL, *jadaPrecRight=NULL;
    
  // the RHS used, if there is no left preconditioning this is just a pointer to res, otherwise
  // it is allocated and the preconditioner is applied
  TYPE(mvec_ptr) rhs = (TYPE(mvec_ptr))res;

  // wrap the preconditioner so that apply_shifted is called
  if (me->leftPrecon!=NULL)
  {
    // create preconditioner with totalNumSys zero shifts because we need to apply it to the RHS as well and
    // we check internally that sufficient shifts are given.
    PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->leftPrecon,q,Bq,&preconShifts[0],totalNumRHS,&jadaPrecL,me->preconSkewProject,iflag),*iflag);
    rhs=NULL;
    PHIST_CHK_IERR(SUBR(mvec_clone_shape)(&rhs,res,iflag),*iflag);
    PHIST_CHK_IERR(jadaPrecL.apply(st::one(),jadaPrecL.A,res,st::zero(),rhs,iflag),*iflag);
    jadaPrecLeft=&jadaPrecL;
    PHIST_CHK_IERR(SUBR(jadaOp_set_leftPrecond)(&jadaOp,jadaPrecLeft,iflag),*iflag);
  }
  else if (me->rightPrecon!=NULL)
  {
    PHIST_CHK_IERR(SUBR(jadaPrec_create)(me->rightPrecon,q,Bq,&preconShifts[0],totalNumSys,&jadaPrecR,me->preconSkewProject,iflag),*iflag);
    jadaPrecRight=&jadaPrecR;
  }
  
  // in case a preconditioned copy was made, delete it at the end of the scope
  MvecOwner<_ST_> _rhs(rhs!=res?rhs:NULL);

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
        if( resIndex != NULL ) ind = resIndex[ind];
        PHIST_CHK_IERR(SUBR(mvec_view_block)(rhs, &res_i, ind, ind, iflag), *iflag);
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
      PHIST_CHK_NEG_IERR(SUBR(blockedMINRESstates_iterate)(&jadaOp, jadaPrecRight, &activeStates[0], k, &nIter, iflag), 
      *iflag);
    }
    else if (me->method_==phist_CARP_CG)
    {
      *iflag=PHIST_NOT_IMPLEMENTED;
    }
    else
    {
      int nIter=maxIter;
      PHIST_CHK_NEG_IERR(SUBR(blockedGMRESstates_iterate)(&jadaOp, jadaPrecRight,&activeStates[0], k, &nIter, useIMGS, iflag), *iflag);
    }

    // optimization to use always full blocks and ignore tolerances
    if( abortAfterFirstConvergedInBlock && k > 0 )
    {
      int ind = index[activeStates[0]->id];
      PHIST_CHK_IERR(SUBR(mvec_view_block)(t, &t_i, ind, ind+k-1, iflag), *iflag);
      _MT_ tmp[k];
      PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(&activeStates[0], k, jadaPrecRight, t_i, tmp, false, iflag), 
      *iflag);
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
          PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(&activeStates[i],1,jadaPrecRight, t_i, &tmp, false, iflag), 
          *iflag);

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

  // delete the jadaOp and preconditioner wrappers
  if (jadaPrecLeft!=NULL)
  {
    PHIST_CHK_IERR(jadaPrecLeft->destroy(jadaPrecLeft, iflag), *iflag);
  }
  if (jadaPrecRight!=NULL)
  {
    PHIST_CHK_IERR(jadaPrecRight->destroy(jadaPrecRight, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR(jadaOp_delete)(&jadaOp, iflag), *iflag);
  *iflag = nUnconvergedSystems;
}
