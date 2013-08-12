void SUBR(op_wrap_crsMat)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A, int* ierr)
  {
  *ierr=0;
  op->A_ = A;
  PHIST_CHK_IERR(SUBR(crsMat_get_range_map)(A,&op->range_map_,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(crsMat_get_domain_map)(A,&op->domain_map_,ierr),*ierr);
  op->apply = &SUBR(crsMat_times_mvec);
  return;
  }

