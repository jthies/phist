//! create a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_create)(TYPE(jadaCorrectionSolver_ptr) *me, int pgmresBlockDim, const_map_ptr_t map, int pgmresMaxBase, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  PHIST_CHK_IERR( *ierr = (pgmresBlockDim <= 0) ? -1 : 0, *ierr);

  *me = new TYPE(jadaCorrectionSolver);
  (*me)->gmresBlockDim_ = pgmresBlockDim;
  (*me)->pgmresStates_  = new TYPE(pgmresState_ptr)[pgmresBlockDim];
  PHIST_CHK_IERR(SUBR(pgmresStates_create)((*me)->pgmresStates_, pgmresBlockDim, map, pgmresMaxBase, ierr), *ierr);
}

//! delete a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_delete)(TYPE(jadaCorrectionSolver_ptr) me, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  PHIST_CHK_IERR(SUBR(pgmresStates_delete)(me->pgmresStates_, me->gmresBlockDim_, ierr), *ierr);
  delete[] me->pgmresStates_;
  delete me;
}


//! calculate approximate solutions to given set of jacobi-davidson correction equations
//!
//! arguments:
//! jdCorrSolver    the jadaCorrectionSolver object
//! A_op            matrix A passed to jadaOp_create
//! B_op            matrix B passed to jadaOp_create
//! Qtil            projection vectors V passed to jadaOp_create
//! BQtil           projection vectors BV passed to jadaOp_create
//! sigma           (pos.!) shifts, -sigma[i], i in {1, ..., nvec} is passed to the jadaOp
//! res             JD residuals, e.g. rhs of the correction equations
//! resIndex        if not NULL, specifies permutation of the residual array to avoid unnecessary copying in the jada-algorithm
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! t               returns approximate solution vectors
//! ierr            a value > 0 indicates the number of systems that have not converged to the desired tolerance
void SUBR(jadaCorrectionSolver_run)(TYPE(jadaCorrectionSolver_ptr) me,
                                    TYPE(const_op_ptr)    A_op,     TYPE(const_op_ptr)    B_op, 
                                    TYPE(const_mvec_ptr)  Qtil,     TYPE(const_mvec_ptr)  BQtil,
                                    const _ST_            sigma[],  TYPE(const_mvec_ptr)  res,      const int resIndex[], 
                                    const _MT_            tol[],    int                   maxIter,
                                    TYPE(mvec_ptr)        t,        int *                 ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  PHIST_CHK_IERR(*ierr = (maxIter <= 0) ? -1 : 0, *ierr);

  // set solution vectors to zero to add them up later
  PHIST_CHK_IERR(SUBR(mvec_put_value)(t, st::zero(), ierr), *ierr);

  // make sure all states are resetted
  for(int i = 0; i < me->gmresBlockDim_; i++)
  {
    PHIST_CHK_IERR(SUBR(pgmresState_reset)(me->pgmresStates_[i], NULL, NULL, ierr), *ierr);
  }

  // current and maximal block dimension
  int max_k = me->gmresBlockDim_;
  int k = max_k;

  // total number of systems to solve
  int totalNumSys;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(t, &totalNumSys, ierr), *ierr);

  // index of currently iterated systems in all systems to solve
  std::vector<int> index(max_k);

  // shifts currently in use
  std::vector<_ST_> currShifts(max_k, st::zero());

  // we need a jadaOp
  TYPE(op) jadaOp;
  PHIST_CHK_IERR(SUBR(jadaOp_create)(A_op, B_op, Qtil, BQtil, &currShifts[0], k, &jadaOp, ierr), *ierr);

  // the next system to consider if one converged/failed
  int nextSystem = 0;

  // just used to prevent inifite loops
  int finiteLoopCounter = 0;

  // we need some views
  TYPE(mvec_ptr) res_i = NULL;
  TYPE(mvec_ptr) t_i   = NULL;

  // total number of pgmres-iterations performed, not used
  int nTotalIter = 0;

  // number of unconverged systems to be returned
  int nUnconvergedSystems = 0;

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
      if( me->pgmresStates_[i]->status == -2 )
      {
        // setup the next system waiting to be solved
        int ind = nextSystem;
        if( resIndex != NULL )
          ind = resIndex[ind];
        PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res, &res_i, ind, ind, ierr), *ierr);
        PHIST_CHK_IERR(SUBR(pgmresState_reset)(me->pgmresStates_[i], res_i, NULL, ierr), *ierr);
        me->pgmresStates_[i]->tol = tol[nextSystem];
        index[me->pgmresStates_[i]->id] = nextSystem;

        nextSystem++;
      }
    }

    // gather systems that are waiting to be iterated
    std::vector<TYPE(pgmresState_ptr)> activeStates;
    int firstId = max_k;
    PHIST_SOUT(PHIST_INFO, "Iterating systems:");
    for(int i = 0; i < max_k; i++)
    {
      if( std::abs(me->pgmresStates_[i]->status) == 1 )
      {
        PHIST_SOUT(PHIST_INFO, "\t%d", index[me->pgmresStates_[i]->id]);
        activeStates.push_back(me->pgmresStates_[i]);
        firstId = std::min(firstId,activeStates.back()->id);
      }
    }
    PHIST_SOUT(PHIST_INFO, "\n");
    k = activeStates.size();

    // set correct shifts
    for(int i = 0; i < k; i++)
    {
      currShifts[activeStates[i]->id - firstId] = -sigma[index[activeStates[i]->id]];
    }

    // actually iterate
    PHIST_CHK_NEG_IERR(SUBR(pgmresStates_iterate)(&jadaOp, &activeStates[0], k, &nTotalIter, ierr), *ierr);


    // check the status of the systems
    for(int i = 0; i < k; i++)
    {
      if( activeStates[i]->status == 0 || activeStates[i]->status == 2 )
      {
        // update solution
        int ind = index[activeStates[i]->id];
        PHIST_CHK_IERR(SUBR(mvec_view_block)(t, &t_i, ind, ind, ierr), *ierr);
        _MT_ tmp;
        PHIST_CHK_IERR(SUBR(pgmresStates_updateSol)(&activeStates[i], 1, t_i, &tmp, false, ierr), *ierr);

        if( activeStates[i]->status == 2 && activeStates[i]->totalIter >= maxIter )
          nUnconvergedSystems++;
        else if( activeStates[i]->status == 2 )
        {
          // prepare restart
          PHIST_CHK_IERR(SUBR(pgmresState_reset)(activeStates[i], NULL, t_i, ierr), *ierr);
          continue;
        }
#ifdef TESTING
{
  // determine real residual for comparison
  currShifts[0] = -sigma[ind];
  TYPE(mvec_ptr) res_i = NULL;
  int resInd = ind;
  if( resIndex != NULL )
    resInd = resIndex[resInd];
  PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))res, &res_i, resInd, resInd, ierr), *ierr);
  PHIST_CHK_IERR(jadaOp.apply(-st::one(), jadaOp.A, t_i, st::one(), res_i, ierr), *ierr);
  _MT_ nrm;
  PHIST_CHK_IERR(SUBR(mvec_norm2)(res_i, &nrm, ierr), *ierr);
  PHIST_SOUT(PHIST_INFO,"est. / exp. residual of system %d: %8.4e / %8.4e\n", ind, tmp, nrm);
}
#endif

        // reset to be free in the next iteration
        PHIST_CHK_IERR(SUBR(pgmresState_reset)(activeStates[i], NULL, NULL, ierr), *ierr);
      }
    }

  }


  // normalize result vectors, TODO: should be done in updateSol/pgmres?
  _MT_ tmp[totalNumSys];
  PHIST_CHK_IERR(SUBR(mvec_normalize)(t, tmp, ierr), *ierr);

  // delete views
  PHIST_CHK_IERR(SUBR(mvec_delete)(res_i, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(t_i,   ierr), *ierr);
  // delete the jadaOp
  PHIST_CHK_IERR(SUBR(jadaOp_delete)(&jadaOp, ierr), *ierr);
}

