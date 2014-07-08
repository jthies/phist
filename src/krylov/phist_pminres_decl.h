//!
//! similar to pgmresStates_iterate, but uses MINRES (e.g. the Lanczos procedure instead of Arnoldi)
//! suited only for symmetric/hermitian matrices
//!
void SUBR( pminresStates_iterate ) (TYPE(const_op_ptr) Op, TYPE(pgmresState_ptr) S_array[], int numSys, int* nIter, int* ierr);
