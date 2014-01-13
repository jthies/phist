// this function can be used to create an operator which encapsulates a CRS matrix.
// It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(op_wrap_crsMat)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A, int* ierr)
  {
  *ierr=0;
  op->A = A;
  PHIST_CHK_IERR(SUBR(crsMat_get_range_map)(A,&op->range_map,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(crsMat_get_domain_map)(A,&op->domain_map,ierr),*ierr);
  op->apply = &SUBR(crsMat_times_mvec);
  return;
  }


//
void SUBR(private_apply_identity)(_ST_ alpha, const void* A, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
        {
        *ierr=0;
        TOUCH(A)
        SUBR(mvec_add_mvec)(alpha,X,beta,Y,ierr);
        }

// setup identity operator that returns Y=alpha*X + beta*Y
void SUBR(op_identity)(TYPE(op_ptr) op, int* ierr)
  {
  *ierr=0;
  op->A=NULL;
  op->apply = &SUBR(private_apply_identity);
  }

