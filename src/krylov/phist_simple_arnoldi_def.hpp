//! a simple Arnoldi process to start up the JaDa iteration.
//! Given a minimum basis size m, compute V(:,1:m+1), 
//! H(1:m+1,1:m) such that A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
//! input: v0, V and H allocated with m+1 resp. m columns
//! and nloc resp. m+1 rows.
//!
//! TODO - block Arnoldi, cf. matlab/krylov/arnoldi.m for a prototype.
//! TODO - we may want to include a check if any Ritz values have already
//!        converged in Arnoldi, but this requires some additional programming
//!        effort which we leave for later.
void SUBR(simple_arnoldi)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, int m, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!\n");
    PHIST_CHK_IERR(*ierr=PHIST_NOT_IMPLEMENTED, *ierr);
  }
  // check dimensions
  {
    int nrH,ncH,ncV;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&ncV,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(H,&nrH,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(H,&ncH,ierr),*ierr);
    // note: we actually don't care about the dimensions of V and H as long as they are large 
    // enough since we just work with views in this function. However, we issue a warning 
    // (positive ierr) if the array dimensions are larger than expected so that the user 
    // doesn't unexpectedly end up with empty rows/cols
    if (ncV<m+1 || nrH<m+1 || ncH<m) {*ierr = -1; return;}
    if (ncV!=m+1 || nrH!=m+1 || ncH!=m)
    {
      PHIST_SOUT(PHIST_VERBOSE,"REMARK: input vectors/matrix to arnoldi are larger than necessary.\n");
      PHIST_SOUT(PHIST_VERBOSE,"        Requested subspace dimension is %d,\n",m);
      PHIST_SOUT(PHIST_VERBOSE,"        H is %dx%d (expecting %dx%d), V with %d cols (expecting %d)\n",
                                        nrH, ncH, m+1,m,ncV,m+1);
    }
  }

  // allocate temporary storage
  TYPE(mvec_ptr) v = NULL, av = NULL;
  PHIST_CHK_IERR(SUBR(mvec_create) (&v,  A_op->domain_map, 1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_create) (&av, A_op->domain_map, 1, ierr), *ierr);

  // views in H and V
  TYPE(sdMat_ptr) R1 = NULL, R2 = NULL;
  TYPE(mvec_ptr) Vprev = NULL;

  // normalize v0
  PHIST_CHK_IERR(SUBR(mvec_set_block) (v, v0, 0, 0, ierr), *ierr);
  _MT_ v0norm;
  PHIST_CHK_IERR(SUBR(mvec_normalize) (v, &v0norm,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_set_block) (V, v,  0, 0, ierr), *ierr);

  // initialize H
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(H,st::zero(),ierr),*ierr);


  // Arnoldi loop
  for(int i = 0; i < m; i++)
  {
    // apply A
    PHIST_CHK_IERR(A_op->apply(st::one(), A_op->A, v, st::zero(), av, ierr), *ierr);
    // copy it to v to improve the accuracy of R1 later
    PHIST_CHK_IERR(SUBR(mvec_set_block) (v,  av, 0, 0, ierr), *ierr);

    // store AV directly if needed later
    if( AV != NULL )
    {
      PHIST_CHK_IERR(SUBR(mvec_set_block) (AV, av, i, i, ierr), *ierr);
    }

    // orthogonalize: Q*R1 = W-VR2
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vprev,0,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R2,0,i,i,i,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,i+1,i+1,i,i,ierr),*ierr);
    PHIST_CHK_NEG_IERR(SUBR(orthog)(Vprev,av,R1,R2,3,ierr),*ierr);
#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_DEB("Arnoldi-step %d\n",i);
    PHIST_CHK_IERR(SUBR(sdMat_print)(H,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_print)(R2,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_print)(R1,ierr),*ierr);
#endif    
    if( *ierr > 0 )
    {
      PHIST_SOUT(PHIST_INFO,"found invariant subspace in arnoldi, expanding basis with a randomly generated orthogonal vector\n");
    }
    else
    {
      // the result for R1 from av'*(A*v) may not be precise enough, so improve it!
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(), av, v, st::zero(), R1, ierr), *ierr);
    }

    // store result in V
    PHIST_CHK_IERR(SUBR(mvec_set_block)(V, av, i+1, i+1, ierr), *ierr);

    // swap vectors
    std::swap(v, av);
  }


  // delete views
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R2, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vprev, ierr), *ierr);

  // delete temp. arrays
  PHIST_CHK_IERR(SUBR(mvec_delete)(v,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(av, ierr), *ierr);
}


//! starts with random blockvector
void SUBR(simple_blockArnoldi)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
                               TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV,
                               TYPE(sdMat_ptr) H, int m, int bs, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;
  if( B_op != NULL )
  {
    PHIST_SOUT(PHIST_ERROR,"case B_op != NULL (e.g. B != I) not implemented yet!\n");
    PHIST_CHK_IERR(*ierr=PHIST_NOT_IMPLEMENTED, *ierr);
  }
  // check dimensions
  {
    int nrH,ncH,ncV;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&ncV,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(H,&nrH,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(H,&ncH,ierr),*ierr);
    // note: we actually don't care about the dimensions of V and H as long as they are large 
    // enough since we just work with views in this function. However, we issue a warning 
    // (positive ierr) if the array dimensions are larger than expected so that the user 
    // doesn't unexpectedly end up with empty rows/cols
    PHIST_CHK_IERR( *ierr =  (ncV<m+bs || nrH<m+bs || ncH<m) ? -1 : 0, *ierr);
    if (ncV!=m+bs || nrH!=m+bs || ncH!=m)
    {
      PHIST_SOUT(PHIST_VERBOSE,"REMARK: input vectors/matrix to arnoldi are larger than necessary.\n");
      PHIST_SOUT(PHIST_VERBOSE,"        Requested subspace dimension is %d,\n",m);
      PHIST_SOUT(PHIST_VERBOSE,"        H is %dx%d (expecting %dx%d), V with %d cols (expecting %d)\n",
                                        nrH, ncH, m+1,m,ncV,m+1);
    }
  }
  // check m is a multiple of the blocksize
  PHIST_CHK_IERR( *ierr = (m % bs != 0) ? -1 : 0, *ierr);


  // allocate temporary storage
  TYPE(mvec_ptr) v = NULL, av = NULL;
  PHIST_CHK_IERR(SUBR(mvec_create) (&v,  A_op->domain_map, bs, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_create) (&av, A_op->domain_map, bs, ierr), *ierr);

  // views in H and V
  TYPE(sdMat_ptr) R1 = NULL, R2 = NULL;
  TYPE(mvec_ptr) Vprev = NULL;

  // create random orthogonal block vector
  PHIST_CHK_IERR(SUBR(mvec_random) (v, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,0,bs-1,0,bs-1,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_QR)(v, R1, ierr), *ierr);
  // copy to V
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V, v, 0, bs-1, ierr), *ierr);

  // initialize H
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(H,st::zero(),ierr),*ierr);


  // Arnoldi loop
  for(int i = 0; i < m/bs; i++)
  {
    // apply A
    PHIST_CHK_IERR(A_op->apply(st::one(), A_op->A, v, st::zero(), av, ierr), *ierr);
    // copy it to v to improve the accuracy of R1 later
    PHIST_CHK_IERR(SUBR(mvec_set_block) (v,  av, 0, bs-1, ierr), *ierr);

    // store AV directly if needed later
    if( AV != NULL )
    {
      PHIST_CHK_IERR(SUBR(mvec_set_block) (AV, av, i*bs, (i+1)*bs-1, ierr), *ierr);
    }

    // orthogonalize: Q*R1 = W-VR2
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vprev,0,(i+1)*bs-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R2,0,(i+1)*bs-1,i*bs,(i+1)*bs-1,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(H,&R1,(i+1)*bs,(i+2)*bs-1,i*bs,(i+1)*bs-1,ierr),*ierr);
    PHIST_CHK_NEG_IERR(SUBR(orthog)(Vprev,av,R1,R2,3,ierr),*ierr);
    if( *ierr > 0 )
    {
      PHIST_SOUT(PHIST_INFO,"found invariant subspace in arnoldi, expanding basis with a randomly generated orthogonal vector\n");
    }
    else
    {
      // the result for R1 from av'*(A*v) may not be precise enough, so improve it!
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(), av, v, st::zero(), R1, ierr), *ierr);
    }

    // store result in V
    PHIST_CHK_IERR(SUBR(mvec_set_block)(V, av, (i+1)*bs, (i+2)*bs-1, ierr), *ierr);

    // swap vectors
    std::swap(v, av);
  }


  // delete views
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R1, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(R2, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vprev, ierr), *ierr);

  // delete temp. arrays
  PHIST_CHK_IERR(SUBR(mvec_delete)(v,  ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_delete)(av, ierr), *ierr);
}
