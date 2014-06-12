//! create a feastCorrectionSolver object
void SUBR(feastCorrectionSolver_create)(TYPE(feastCorrectionSolver_ptr) *me, 
        TYPE(const_crsMat_ptr) A,
        linSolv_t method,
        int blockSize,
        int numShifts, _MT_ sigma_r[], _MT_ sigma_i[],
        _MT_ tol, int maxIter,
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
    const_map_ptr_t map=NULL;
    PHIST_CHK_IERR(SUBR(crsMat_get_row_map)(A,&map,ierr),*ierr);
    // create one CARP-CG object per shift.
    (*me)->carp_cgStates_ = new TYPE(carp_cgState_ptr)[numShifts];
    
    PHIST_CHK_IERR(SUBR(carp_cgStates_create)
        ((*me)->carp_cgStates_,numShifts,map,blockSize,maxIter,ierr),*ierr);
    
    for (int i=0;i<numShifts;i++)
    {
      (*me)->carp_cgStates_[i]->tol = tol;
      (*me)->carp_cgStates_[i]->sigma_r = sigma_r[i];
      (*me)->carp_cgStates_[i]->sigma_i = sigma_i[i];
    }
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
//! fCorrSolver    the feastCorrectionSolver object
//! rhs            rhs of the correction equations, rhs is the same for all
//!                to shift[i] in create()
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! sol             returns approximate solution vectors, sol[i] belongs to shift[i] and rhs
//! ierr            a value > 0 indicates the number of systems that have not converged to the 
//!                 desired tolerance
void SUBR(feastCorrectionSolver_run)(TYPE(feastCorrectionSolver_ptr) fCorrSolver,
                                    TYPE(const_mvec_ptr) rhs,
                                    TYPE(mvec_ptr) sol_r[],
                                    TYPE(mvec_ptr) sol_i[],
                                    int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (fCorrSolver->method_==CARP_CG)
  {
    *ierr=-99;
  }
  else
  {
    *ierr=-99;
  }
}
