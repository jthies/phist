
// these are some thoughts I wrote down some time ago, I keep
// them here for the moment to look at later

//! in Jacobi-Davidson, we want ot perform an iterative solve
//! (GMRES to start with) for the correction equation. The   
//! corresponding operator we put in the GMRES solver is     
//!                                                          
//! (I-BVV') M\(A-sigma_j B),                                 
//!                                                          
//! where V = [Q v] is the basis                             
//! of converged eigenvectors combined with the current      
//! search space, and M is some preconditioner for a nearby  
//! matrix A-tau*B). We want the GMRES method to run in para-
//! lel for several shifts sigma_j, which we achieve by using
//! the task-buffer programming model (cf. examples/sched/...
//! main_task_model.c for a simple example). For now we let  
//! the jadaOp handle the orthogonalization wrt. V, but in   
//! the end we can hopefully use some bordered preconditio-  
//! ning (like in HYMLS) for the matrix                      
//!                                                          
//! A-tau*B  BV                                               
//!     V'   0                                               
//!                                                          
//! There are in principle two options, use (A-sigma_jI) as  
//! operator and "(I-BVV')M\" as left preconditioner, or use  
//! (I-BVV')(A-sigma_jB) as operator and M\ as right precond. 
//! The first is more suitable for the case where the pre-   
//! conditioner can handle the border, the second if flexible
//! ible GMRES has to be used, e.g. in case of nested ite-   
//! rations on the separators in DSC. Combining the right    
//! preconditioner with the projection should also be fine,  
//! the operator is then (A-sigma_j B) M\(I-VV'B), I think    
//! that should give the same results more or less.
//!
//! Remark: for the generalized eigenvalue problem with hpd. B
//!         we need always both projections, e.g. Y = (I-BVV')(A*X_-B*X_*sigma)
//!         with X_ = (I-VV'B)X
//!         (melven)
//! Remark: if we directly store -sigma, we can avoid an additional vector operation
//!         (as the mvec_vadd_mvec doesn't provei an additional scaling argument)


// private struct to keep all the pointers we need in order to apply the operator.
typedef struct TYPE(jadaOp_data)
{
  TYPE(const_op_ptr)    A_op;   // operator of the general matrix A
  TYPE(const_op_ptr)    B_op;   // operator of the hpd. matrix B, assumed I when NULL
  TYPE(const_mvec_ptr)  V;      // B-orthonormal basis
  TYPE(const_mvec_ptr)  BV;     // B*V
  const _ST_*           sigma;  // array of NEGATIVE shifts, assumed to have correct size; TODO: what about 'complex' shifts for real JDQR?
  TYPE(sdMat_ptr)       VY;     // temporary storage for V' times intermediate Y for the projection (I-VV')Y
  TYPE(mvec_ptr)        AX;     // temporary storage for A*X, respectively  A*X_proj for B!=NULL
  TYPE(mvec_ptr)        BX;     // temporary storage for B*X_proj, only needed for B != NULL
  TYPE(mvec_ptr)        X_proj; // temporary storage for (I-VV'B)X, only used for B!= NULL
} TYPE(jadaOp_data);


// actually applies the jada-operator:
// for B==NULL:  Y <- alpha* (I-VV')(AX+BX*sigma) + beta*Y
// for B!=NULL:  Y <- alpha* (I-BVV')*(A(I-VV'B)X + BX*sigma) + beta*Y
void SUBR(jadaOp_apply)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *ierr);

  TYPE(const_mvec_ptr) BX;
  if( jadaOp->B_op == NULL )
  {
    // calculate AX and set BX = 
    PHIST_CHK_IERR( jadaOp->A_op->apply(ONE, jadaOp->A_op->A, X, ZERO, jadaOp->AX, ierr), *ierr);               // AX     <- A*X
    BX = X;
  }
  else // B_op != NULL
  {
    // calculate X_proj, AX_proj and BX_proj
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (ONE, jadaOp->BV, X, ZERO, jadaOp->VY, ierr), *ierr);              // VY     <- (BV)'X
    PHIST_CHK_IERR( SUBR( mvec_add_mvec    ) (ONE, X, ZERO, jadaOp->X_proj, ierr), *ierr);                      // X_proj <- X
    PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-ONE, jadaOp->V, jadaOp->VY, ONE, jadaOp->X_proj, ierr), *ierr);  // X_proj <- X_proj - V*VY
    PHIST_CHK_IERR( jadaOp->A_op->apply(ONE, jadaOp->A_op->A, jadaOp->X_proj, ZERO, jadaOp->AX, ierr), *ierr);  // AX     <- A*X_proj
    PHIST_CHK_IERR( jadaOp->A_op->apply(ONE, jadaOp->B_op->A, jadaOp->X_proj, ZERO, jadaOp->BX, ierr), *ierr);  // BX     <- B*X_proj
    BX = jadaOp->BX;
  }


  // Y <- alpha* (I-BVV')(AX-BX*sigma) + beta*X
  PHIST_CHK_IERR( SUBR( mvec_add_mvec    ) (ONE, jadaOp->AX, ZERO, Y, ierr), *ierr);                          // Y      <- AX
  PHIST_CHK_IERR( SUBR( mvec_vadd_mvec   ) (jadaOp->sigma, BX, ONE, Y, ierr), *ierr);                         // Y      <- Y + BX*sigma
  PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (ONE, jadaOp->V, Y, ZERO, jadaOp->VY, ierr), *ierr);               // VY     <- V'*Y
  PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-ONE, jadaOp->BV, jadaOp->VY, ONE, Y, ierr), *ierr);              // Y      <- Y - BV*VY

  if( alpha != ONE || beta != ZERO )
  {
    PHIST_CHK_IERR( SUBR( mvec_add_mvec  ) (beta, X, alpha, Y, ierr), *ierr);
  }
}


// allocate and initialize the jadaOp struct
void SUBR(jadaOp_create)(TYPE(const_op_ptr)    A_op,    TYPE(const_op_ptr)    B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], TYPE(sdMat_ptr)       VY,
                         TYPE(mvec_ptr)        AX,      TYPE(mvec_ptr)        BX,
                         TYPE(mvec_ptr)        X_proj,  TYPE(op_ptr)          jdOp,
                         int*                  ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  // allocate jadaOp struct
  TYPE(jadaOp_data) *myOp = (TYPE(jadaOp_data)*)malloc(sizeof(TYPE(jadaOp_data)));

  // setup jadaOp members
  myOp->A_op   = A_op;
  myOp->B_op   = B_op;
  myOp->V      = V;
  myOp->BV     = (B_op != NULL ? BV     : V);
  myOp->sigma  = sigma;
  myOp->VY     = VY;
  myOp->AX     = AX;
  myOp->BX     = (B_op != NULL ? BX     : NULL);
  myOp->X_proj = (B_op != NULL ? X_proj : NULL);

  // setup op_ptr function pointers
  jdOp->A     = (const void*)myOp;
  jdOp->apply = (&SUBR(jadaOp_apply));
}


// deallocate jadaOp struct
void SUBR(jadaOp_delete)(TYPE(op_ptr) jdOp, int *ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  // get jadaOp
  TYPE(jadaOp_data) *jadaOp = (TYPE(jadaOp_data)*) jdOp->A;

  // delete jadaOp
  free(jadaOp);
}

//! access AX from last call to apply (return a view to it)
void SUBR(jadaOp_view_AX)(TYPE(const_op_ptr) op, TYPE(mvec_ptr)*AX, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *ierr);

  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(jadaOp->AX, &nvec, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(jadaOp->AX, AX, 0, nvec-1, ierr), *ierr);
}

//! access BX from last call to apply (return a view to it)
void SUBR(jadaOp_view_BX)(TYPE(const_op_ptr) op, TYPE(mvec_ptr)*BX, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *ierr);

  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(jadaOp->BX, &nvec, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(jadaOp->BX, BX, 0, nvec-1, ierr), *ierr);
}

//! access X_proj from last call to apply (return a view to it)
void SUBR(jadaOp_view_X_proj)(TYPE(const_op_ptr) op, TYPE(mvec_ptr)*X_proj, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *ierr);

  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(jadaOp->X_proj, &nvec, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_view_block)(jadaOp->X_proj, X_proj, 0, nvec-1, ierr), *ierr);
}

