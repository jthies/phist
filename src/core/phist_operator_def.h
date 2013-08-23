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

//! this is a helper function. Given an array of multi-vectors, 
//! calculates the displacement between the first element of each
//! column of [v[0],v[1],v[2],...]). If the displacement is not
//! constant, returns -1. For contiguously stored multi-vectors,
//! returns the LDA, the leading dimension of the data arrays, which
//! are in that case to be all the same.
void SUBR(check_stride)(TYPE(const_mvec_ptr) v[], int n_mvec,
        int* stride, int* nvec_tot, int* ierr)
  {
  int i;
  int nvec_i;
  lidx_t lda_i;
  _ST_ *val_i, *val_im1;
  *ierr=0;
  *stride=0;
  *nvec_tot=0;
  if (n_mvec==0) return;
  
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v[0],&nvec_i,ierr),*ierr);
  *nvec_tot=nvec_i;
  PHIST_CHK_IERR(SUBR(mvec_extract_view)(v[0],&val_im1,stride,ierr),*ierr);

  for (i=1; i<n_mvec;i++)
    {
    if (*stride>0)
      {
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(v[i],&val_i,&lda_i,ierr),*ierr);
      if (lda_i!=*stride)
        {
        *stride=-1;
        }
      if (val_i!=val_im1+nvec_i*lda_i)
        {
        *stride=-2;
        }
      }
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v[i],&nvec_i,ierr),*ierr);
    *nvec_tot+=nvec_i;
    }
  return;
  }

// apply operator to several multi-vectors, which are not necessarily contiguous in memory.
void SUBR(op_apply_buffered)(argList_t* args)
  {
  int i;
  int n_mvecs=args->n;
  int* ierr=&args->ierr;
  TYPE(const_op_ptr) A_op = (TYPE(const_op_ptr))args->shared_arg;
  TYPE(const_mvec_ptr)* X = (TYPE(mvec_ptr)*)args->in_arg;
  TYPE(mvec_ptr)* Y = (TYPE(mvec_ptr)*)args->out_arg;
  
  TYPE(mvec_ptr) myX;
  TYPE(mvec_ptr) myY;

  const_map_ptr_t map;
  int contig;
  int stride_X, stride_Y;
  int nvec_tot_X, nvec_tot_Y;
  int nvec_X, nvec_Y;
  _ST_ *val_X, *val_Y;

  int imin,imax;

  *ierr=0;

  // no input?
  if (n_mvecs==0) return;
  
  contig=1; // indicates if vectors lie in memory with constant stride
  
  // just one mvec?
  if (n_mvecs==1)
    {
    myX=X[0];
    myY=Y[0];
    }
  else
    {
    PHIST_CHK_IERR(SUBR(mvec_get_map)(X[0],&map,ierr),*ierr);

    // count vectors and check if vectors lie in memory with 
    // constant stride => create a view (if contig=1 after this loop)
    // TODO - once we have established this we should skip the test
    // until the setup changes.
      PHIST_CHK_IERR(SUBR(check_stride)(X,n_mvecs,&stride_X,&nvec_tot_X,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(check_stride)(Y,n_mvecs,&stride_Y,&nvec_tot_Y,ierr),*ierr);

    // inconsistent input - not the same number of vectors in total
    PHIST_CHK_IERR(*ierr=-(nvec_tot_X==nvec_tot_Y),*ierr);
    contig = (stride_X>0 && stride_X==stride_Y)? 1:0;

    if (contig)
      {
      // vectors lie in memory with constant stride => can
      // create a view and use e.g. BLAS 3 operations on the whole block
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(X[0],&val_X,&stride_X,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_extract_view)(Y[0],&val_Y,&stride_Y,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_create_view)(&myX,map,val_X,stride_X,nvec_tot_X,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_create_view)(&myY,map,val_Y,stride_Y,nvec_tot_Y,ierr),*ierr);
      }
    else
      {
      // create a copy and copy vectors into adjacent memory
      PHIST_CHK_IERR(SUBR(mvec_create)(myX,map,nvec_tot_X,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_create)(myY,map,nvec_tot_Y,ierr),*ierr);
      imin=0;
      for (i=0;i<n_mvecs;i++)
        {
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X[i],&nvec_X,ierr),*ierr);
        PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Y[i],&nvec_Y,ierr),*ierr);
        // check that the vectors X_i and Y_i have the same number of columns,
        // this hasn't been checked in check_stride, but it makes sense to    
        // expect that, even though it is not strictly necessary.
        PHIST_CHK_IERR(*ierr=(nvec_X!=nvec_Y),*ierr);
        imax=imin+nvec_X-1;
        PHIST_CHK_IERR(SUBR(mvec_set_block)(myX,X[i],imin,imax,ierr),*ierr);
        imin=imax+1;
        }
      }
    }
        
  // perform the operation
  PHIST_CHK_IERR(A_op->apply(ONE,A_op,myX,ZERO,myY,ierr),*ierr);
  
  if (!contig)
    {
    // copy the vectors back to where they belong
    imin=0;    
    for (i=0;i<n_mvecs;i++)
      {
      PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X[i],&nvec_X,ierr),*ierr);
      imax=imin+nvec_X-1;
      PHIST_CHK_IERR(SUBR(mvec_get_block)(myY,Y[i],imin,imax,ierr),*ierr);
      imin=imax+1;
      }
    
    PHIST_CHK_IERR(SUBR(mvec_delete)(myX,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_delete)(myY,ierr),*ierr);
    }// not contiguous, we have to copy data
  
  return;
  }

//
void SUBR(private_apply_identity)(_ST_ alpha, TYPE(const_op_ptr) A, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
        {
        SUBR(mvec_add_mvec)(alpha,X,beta,Y,ierr);
        }

// setup identity operator that returns Y=alpha*X + beta*Y
void SUBR(op_identity)(TYPE(op_ptr) op, int* ierr)
  {
  op->A=NULL;
  op->apply = &SUBR(private_apply_identity);
  }

