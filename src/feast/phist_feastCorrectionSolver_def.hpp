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
  
  PHIST_SOUT(PHIST_VERBOSE,"creating feastCorrection solver with\n"
                           " block size %d\n and %d shifts\n",blockSize,numShifts);
  
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
#ifdef PHIST_KERNEL_LIB_GHOST
// GHOST has mixed arithmetic built in, we could simply use the complex kernels
// (with complex shifts and vectors) but pass in a real-valued matrix (in theory...)
#warning "for ghost we should use the complex carp_cg but pass in the real matrix."\
         "Using the real arithmetic version may cost performance..."

#endif  
    // create one CARP-CG object per shift.
    (*me)->carp_cgStates_ = new TYPE(carp_cgState_ptr)[numShifts];
    
    PHIST_CHK_IERR(SUBR(carp_cgStates_create)
        ((*me)->carp_cgStates_,numShifts,(*me)->sigma_r_,(*me)->sigma_i_,
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
  
  // this function solves nshifts*nrhs linear systems and groups them
  // somehow into blocks of size bs. The outer loop if over the columns
  // (of sol and rhs), the inner over the shifts.
  
  int blockSize=me->blockSize_;
  int nrhs;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(rhs,&nrhs,ierr),*ierr);

  bool copy_input_vecs= (nrhs!=blockSize);

  for (int c0=0;c0<nrhs;c0+=blockSize)
  {
  int c1=std::min(c0+blockSize,nrhs)-1;
  int bs=c1-c0+1;
  PHIST_SOUT(PHIST_VERBOSE,"SOLVE SYSTEMS (%d:%d)\n",c0,c1);
    if (me->method_==CARP_CG)
    {
    TYPE(mvec_ptr) b,*x_r,*x_i;
    b=(TYPE(mvec_ptr))rhs;
    x_r=sol_r;
    x_i=sol_i;
      
    if (copy_input_vecs)
    {
      PHIST_SOUT(PHIST_DEBUG,"create tmp vectors\n");
      const_map_ptr_t map;
      PHIST_CHK_IERR(SUBR(mvec_get_map)(b,&map,ierr),*ierr);

      x_r=new TYPE(mvec_ptr)[me->numShifts_];
      x_i=new TYPE(mvec_ptr)[me->numShifts_];

      PHIST_CHK_IERR(SUBR(mvec_create)(&b,map,bs,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_get_block)(rhs,b,c0,c1,ierr),*ierr);
      for (int i=0;i<me->numShifts_;i++)
      {
        PHIST_CHK_IERR(SUBR(mvec_create)(&x_r[i],map,bs,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_create)(&x_i[i],map,bs,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_get_block)(sol_r[i],x_r[i],c0,c1,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_get_block)(sol_i[i],x_i[i],c0,c1,ierr),*ierr);
      }
    }
  
      MT* normsB=NULL;
      // adjust block size if fewer systems are solved.
      if (bs!=blockSize)
      {
        PHIST_SOUT(PHIST_DEBUG,"re-construct CARP structs\n");
        PHIST_CHK_IERR(SUBR(carp_cgStates_delete)
            (me->carp_cgStates_,me->numShifts_,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(carp_cgStates_create)
            (me->carp_cgStates_,me->numShifts_,me->sigma_r_,me->sigma_i_,
            me->A_, bs,ierr),*ierr);    
      }
      // reset all CG states. Use the given sol vectors as starting guess.
      for (int s=0; s<me->numShifts_; s++)
      {
        PHIST_CHK_IERR(SUBR(carp_cgState_reset)(me->carp_cgStates_[s],b,normsB,ierr),*ierr);
        // compute ||b|| only for first system, then copy it
        normsB=me->carp_cgStates_[0]->normB_;
      }
      // now iterate the systems with the different shifts
      // (and the same multiple RHS each). At this point
      // we put the further parallelisation into the hand of
      // the carp_cg implementation, which might delegate shifts
      // to other nodes, use a queuing system etc. The starting
      // vectors are taken as received by this function and passed
      // on to carp_cg.
      PHIST_CHK_IERR(SUBR(carp_cgStates_iterate)
                (me->carp_cgStates_, me->numShifts_, 
                x_r, x_i, tol, maxIter,ierr),*ierr);
      // reset to original block size if fewer systems were solved.
      if (bs!=blockSize)
      {
        PHIST_SOUT(PHIST_DEBUG,"re-construct CARP structs with original block size\n");
        PHIST_CHK_IERR(SUBR(carp_cgStates_delete)
            (me->carp_cgStates_,me->numShifts_,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(carp_cgStates_create)
            (me->carp_cgStates_,me->numShifts_,me->sigma_r_,me->sigma_i_,
            me->A_, blockSize,ierr),*ierr);    
      }
      if (copy_input_vecs)
      {
        PHIST_SOUT(PHIST_DEBUG,"delete tmp vectors\n");
        PHIST_CHK_IERR(SUBR(mvec_delete)(b,ierr),*ierr);
        for (int i=0;i<me->numShifts_;i++)
        {
          PHIST_CHK_IERR(SUBR(mvec_set_block)(sol_r[i],x_r[i],c0,c1,ierr),*ierr);
          PHIST_CHK_IERR(SUBR(mvec_set_block)(sol_i[i],x_i[i],c0,c1,ierr),*ierr);
          PHIST_CHK_IERR(SUBR(mvec_delete)(x_r[i],ierr),*ierr);
          PHIST_CHK_IERR(SUBR(mvec_delete)(x_i[i],ierr),*ierr);
        }
        delete [] x_r;
        delete [] x_i;
      }
    }
    else
    {
      *ierr=-99;
    }
  }
}
