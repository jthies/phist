extern "C" {

// this function can be used to create an operator which encapsulates a CRS matrix.
// It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(op_wrap_sparseMat)(TYPE(op_ptr) op, TYPE(const_sparseMat_ptr) A, int* iflag)
  {
  *iflag=0;
  op->A = A;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&op->range_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&op->domain_map,iflag),*iflag);
  op->apply = &SUBR(sparseMat_times_mvec);
  op->applyT = &SUBR(sparseMatT_times_mvec);
  op->apply_shifted = &SUBR(sparseMat_times_mvec_vadd_mvec);
  return;
  }


//
void SUBR(private_idOp_apply)(_ST_ alpha, const void* A, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
  *iflag=0;
  PHIST_TOUCH(A)
  SUBR(mvec_add_mvec)(alpha,X,beta,Y,iflag);
}

//
void SUBR(private_idOp_apply_shifted)(_ST_ alpha, const void* A, _ST_ const *sigma, 
        TYPE(const_mvec_ptr) X,  _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
  *iflag=0;
  PHIST_TOUCH(A)
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,iflag),*iflag);
  _ST_ shifts[nvec];
  for (int i=0;i<nvec;i++) shifts[i]=sigma[i]*alpha;
  SUBR(mvec_vadd_mvec)(shifts,X,beta,Y,iflag);
}

// setup identity operator that returns Y=alpha*X + beta*Y
void SUBR(op_identity)(TYPE(op_ptr) op, int* iflag)
{
  *iflag=0;
  op->A=NULL;
  op->apply = &SUBR(private_idOp_apply);
  op->applyT = &SUBR(private_idOp_apply);
  op->apply_shifted = &SUBR(private_idOp_apply_shifted);
}

} // extern "C"
