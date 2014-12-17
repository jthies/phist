// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182).
// If the input vector has m columns and rank r, ierr=m-r is
// returned. Columns 0:ierr-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. The third argument, E, 
// should be preallocated by the user with m elements, m being
// the number of columns in V. On successful return (ierr>=0),
// e[j] indicates the norm of V(:,j) before the orthogonali-  
// zation step. On exit, V(:,j), j>*ierr has 2-norm 1.
void SUBR(svqb)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, _MT_* e, int* ierr);
