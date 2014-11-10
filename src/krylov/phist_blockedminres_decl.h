//! \addtogroup linear_solvers
//@{

//! \defgroup minres MinRes solver for symmetric indefinite matrices
//@{

//! perform minres iterations (for symmetric indefinite matrices)

//! similar to blockedGMRESstates_iterate, but uses MINRES (e.g. the Lanczos procedure instead of Arnoldi)
//! suited only for symmetric/hermitian matrices
void SUBR( blockedMINRESstates_iterate ) (TYPE(const_op_ptr) Op, TYPE(blockedGMRESstate_ptr) S_array[], int numSys, int* nIter, int* ierr);

//@}

//@}
