//! create a feastCorrectionSolver object
void SUBR(feastCorrectionSolver_create)(TYPE(feastCorrectionSolver_ptr) *me, 
        TYPE(const_crsMat_ptr) A,
        linSolv_t method,
        int blockSize,
        int numShifts, _MT_ sigma_r[], _MT_ sigma_i[],
        int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
  
  *me = new TYPE(feastCorrectionSolver);
  (*me)->A_=A;
  (*me)->numShifts_=numShifts;
  (*me)->blockSize_=blockSize;
  (*me)->sigma_r_ = new MT[numShifts];
  (*me)->sigma_i_ = new MT[numShifts];
  (*me)->rhs_=NULL;
  (*me)->method_=method;

  for (int i=0;i<numShifts;i++)
  {
    (*me)->sigma_r_[i]=sigma_r[i];
    (*me)->sigma_i_[i]=sigma_i[i];
  }
  
  if (method==CARP_CG)
  {
    // create one CARP-CG object per shift.
    (*me)->carp_cgStates_ = new TYPE(carp_cgState_ptr)[numShifts];
    
    PHIST_CHK_IERR(SUBR(carp_cgStates_create)
        ((*me)->carp_cgStates_,numShifts,sigma_r,sigma_i,
        A, blockSize,ierr),*ierr);    
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "method %d (%s) not implemented",(int)method, linSolv2str(method));
    *ierr=-99;
  }
}

//! delete a feastCorrectionSolver object
void SUBR(feastCorrectionSolver_delete)(TYPE(feastCorrectionSolver_ptr) me, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  me->A_=NULL;
  delete [] me->sigma_r_;
  delete [] me->sigma_i_;
  me->rhs_=NULL;
  if (me->method_==CARP_CG)
  {
    PHIST_CHK_IERR(SUBR(carp_cgStates_delete)(me->carp_cgStates_,me->numShifts_,ierr),*ierr);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "method %d (%s) not implemented",(int)(me->method_), 
        linSolv2str(me->method_));
    *ierr=-99;
  }

  delete me;
}

//! calculate approximate solutions to given set of FEAST correction equations
//!
//! arguments:
//! me    the feastCorrectionSolver object
//! rhs            rhs of the correction equations, rhs is the same for all
//!                to shift[i] in create()
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! sol             returns approximate solution vectors, sol[i] belongs to shift[i] and rhs
//! ierr            a value > 0 indicates the number of systems that have not converged to the 
//!                 desired tolerance
void SUBR(feastCorrectionSolver_run)(TYPE(feastCorrectionSolver_ptr) me,
                                    TYPE(const_mvec_ptr) rhs,
                                    _MT_ tol, int maxIter,
                                    TYPE(mvec_ptr) sol_r[],
                                    TYPE(mvec_ptr) sol_i[],
                                    int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (me->method_==CARP_CG)
  {
    MT* normsB=NULL;
    // reset all CG states. Use the given sol vectors as starting guess.
    for (int i=0; i<me->numShifts_; i++)
    {
      PHIST_CHK_IERR(SUBR(carp_cgState_reset)(me->carp_cgStates_[i],rhs,normsB,ierr),*ierr);
      // compute ||b|| only for first system, then copy it
      normsB=me->carp_cgStates_[0]->normB_;
    }
    // now iterate the systems with the different shifts
    // (and the same multiple RHS each). At this point
    // we put the further parallelisation into the hand of
    // the carp_cg implementation, which might delegate shifts
    // to other nodes, use a queuing system etc.
    int nIter=0; // not sure why exactly this counter, but we'll see.
    PHIST_CHK_IERR(SUBR(carp_cgStates_iterate)
                (me->carp_cgStates_, me->numShifts_, 
                sol_r, sol_i, tol, maxIter,ierr),*ierr);
  }
  else
  {
    *ierr=-99;
  }
}
