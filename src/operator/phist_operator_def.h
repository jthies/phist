void _SUBR_(op_wrap_crsMat)(_TYPE_(op_ptr) op, _TYPE_(const_crsMat_ptr) A, int* ierr)
  {
  *ierr=0;
  op->A_ = A;
  _PHIST_ERROR_HANDLER_(_SUBR_(crsMat_get_range_map)(A,&op->range_map_,ierr),*ierr);
  _PHIST_ERROR_HANDLER_(_SUBR_(crsMat_get_domain_map)(A,&op->domain_map_,ierr),*ierr);
  op->apply = &_SUBR_(crsMat_times_mvec);
  return;
  }

