// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182).
// If the input vector has m columns and rank r, ierr=m-r is
// returned. Columns 0:ierr-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. It is possible to obtain
// a QR factorization by e.g. computing [q,r]=qr(B), R=r\q, at
// least if V is of full rank.
void SUBR(svqb)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, int* ierr);
