// forward declaration of apply function (does not appear in the header because
// we don't want anyone to use it from the outside)
void SUBR(private_msOp_apply)(_ST_ alpha, TYPE(const_op_ptr) msOp, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr);

//! create an operator that behaves like a multi-shift operator for a given
//! base operator A_op. The shifts are copied and their given order defines
//! the default order in which they are applied.
void SUBR(msOp_create)(TYPE(op_ptr)* msOp, TYPE(const_op_ptr) A_op, 
        int num_shifts, _ST_* shifts,int* ierr)
  {
  int i;
  msOp = (TYPE(op_ptr))malloc(sizeof(TYPE(op_ptr)));
  msOp->A_ = malloc(sizeof(TYPE(TYPE(multishift))));
  TYPE(multishift) *ms = (TYPE(multishift)*)msOp->A_;

  ms->shifts__=(_ST_*)malloc(num_shifts*sizeof(_ST_)); 
  ms->num_shifts_=num_shifts; 
  ms->shift_order_=(int*)malloc(num_shifts*sizeof(int));;
  for (i=0;i<num_shifts;i++)
    {
    ms->shift_[i]=shifts[i];
    }
  (*msOp)->apply = &SUBR(private_msOp_apply);
  PHIST_CHK_IERR(SUBR(msOp_reset_order)(*msOp,ierr),*ierr);
  }
  
//! set entry in the ordering array such that given shift j is applied to vector column i
void SUBR(msOp_set_order)(TYPE(op_ptr) msOp, int i, int j, int* ierr)
  {
  TYPE(multishift) *ms = (TYPE(multishift)*)msOp->A_;
  PHIST_CHK_IERR(*ierr=(i<0 || j<0 || i>ms->num_shifts_ || j>ms_num_shifts_)?-1:0,*ierr);
  ms->shift_order_[i]=j;  
  }

//! reset the order in which the shifts are applied to the default, 0:num_shifts-1
void SUBR(msOp_reset_order)(TYPE(op_ptr) msOp, int* ierr)
  {
  TYPE(multishift) *ms = (TYPE(multishift)*)msOp->A_;
  int i;
  for (i=0;i<ms->num_shifts_;i++)
    {
    ms->shift_order_[i]=i;
    }
  }
  
//! set entry in the ordering array such that given shift j is applied to vector column i
void SUBR(msOp_delete)(TYPE(op_ptr) msOp, int* ierr)
  {
  TYPE(multishift) *ms = (TYPE(multishift)*)msOp->A_;
  free(ms->shifts_);
  free(ms->shift_order);
  free(msOp);
  }

//! apply Y_j = alpha*(A + shift_[so[j]] I)*X_j + beta*Y_j with given shift order 'so'.
//! TODO: this should be done in a single step, maybe we could have 
//! a corresponding function in the op_t struct and use it if not NULL.
void SUBR(private_msOp_apply)(_ST_ alpha, TYPE(const_op_ptr) msOp, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
  {
  int i,nvec;
  TYPE(multishift)* ms = (const TYPE(multishift)*)msOp->A_;
  _ST_* shifts;
  PHIST_CHK_IERR(SUBR(mvec_get_num_vectors)(X,&nvec,ierr),*ierr);
  if (alpha!=ZERO)
    {
    int k = MIN(nvec,ms->num_shifts);
    shifts = (_ST_*)malloc(nvec*sizeof(_ST));
    for (i=0;i<k;i++)
      {
      shifts[i] = alpha*ms->shifts_[ms->shift_order_[i]];
      }
    for (i=k;i<nvec;i++)
      {
      shifts[i]=ZERO;
      }
    // first set Y_j = beta*Y_j + alpha*shift_j X_j
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(shifts,X,beta,Y,ierr),*ierr);
    // then add Y += alpha*A*X
    PHIST_CHK_IERR(ms->A_op_->apply(alpha,ms->A_op,X,ONE,Y);
    free(shifts);
    }
  else
    {
    // alpha=0 => just compute Y=beta*Y
    PHIST_CHK_IERR(SUBR(mvec_scale)(Y,beta,ierr),*ierr);
    }
  }
