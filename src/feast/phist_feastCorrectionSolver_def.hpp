#ifdef IS_DOUBLE
#define MSUBR(sub) phist_D ## sub
#else
#define MSUBR(sub) phist_S ## sub
#endif
//! create a feastCorrectionSolver object
void SUBR(feastCorrectionSolver_create)(TYPE(feastCorrectionSolver_ptr) *me, 
        TYPE(const_crsMat_ptr) A,
        int blockSize,
        int numShifts, _MT_ sigma_r[], _MT_ sigma_i[],
        _MT_ tol, int maxIter,
        int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
  
  if (numShifts!=1)
  {
    PHIST_SOUT(PHIST_ERROR,"we only allow one shift at a time in this first \n"
                          "prototype of the feastCorrectionSolver, numShifts=1.\n"
                          "You have to create a separate object for each shift and\n"
                          "run the solvers separately. (file %s, line %d)",__FILE__,__LINE__);
    *ierr=-99;
    return;
  }
  
  *me = new TYPE(feastCorrectionSolver);
  (*me)->A_=A;
  (*me)->numShifts_=numShifts;
  (*me)->blockSize_=blockSize;
  (*me)->sigma_r_ = new MT[numShifts];
  (*me)->sigma_i_ = new MT[numShifts];
  (*me)->rhs_=NULL;
  for (int i=0;i<numShifts;i++)
  {
    (*me)->sigma_r[i]=sigma_r[i];
    (*me)->sigma_i[i]=sigma_i[i];
  }
  if (method==CARP_CG)
  {
    phist_map_t* map=NULL;
    PHIST_CHK_IERR(SUBR(crsMat_get_row_map)(A,&map,ierr),*ierr);
    // create one CARP-CG object per shift and rhs vector, leave it to the
    // carp-cg implementation how this is solved and parallelized.
    int numSys = numShifts*blockSize;
    PHIST_CHK_IERR(SUBR(carp_cgStates_create)(&((*me)->carp_cgStates_),numSys,
                map, tol, maxIter,ierr),*ierr);            
    for (int i=0;i<numShifts;i++)
    {
      for (int j=0;j<blockSize;j++
      {
        (*me)->carp_cgStates_[j*blockSize+i]->sigma_r = sigma_r[i];
        (*me)->carp_cgStates_[j*blockSize+i]->sigma_i = sigma_i[i];
      }
    }
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "method %d (%s) not implemented",(int)method, linSolv2str(method));
    *ierr=-99;
  }
}

//! delete a jadaCorrectionSolver object
void SUBR(feastCorrectionSolver_delete)(TYPE(feastCorrectionSolver_ptr) me, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  (*me)->A_=NULL;
  delete [] (*me)->sigma_r_;
  delete [] (*me)->sigma_i_;
  (*me)->rhs_=NULL;
  if (method==CARP_CG)
  {
    int numSys = me->numShifts_ * me->blockSize_;
    PHIST_CHK_IERR(SUBR(carp_cgStates_delete)((*me)->carp_cgStates_,ierr),*ierr);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "method %d (%s) not implemented",(int)method, linSolv2str(method));
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

  if (numShifts!=1)
  {
    PHIST_SOUT(PHIST_ERROR,"we only allow one shift at a time in this first \n"
                          "prototype of the feastCorrectionSolver, numShifts=1.\n"
                          "You have to create a separate object for each shift and\n"
                          "run the solvers separately. (file %s, line %d)",__FILE__,__LINE__);
    *ierr=-99;
    return;
  }

  *ierr=-99;
}
