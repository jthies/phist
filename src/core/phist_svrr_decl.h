// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
// with rank-revealing pivoted cholesky (SVRR)
// If the input vector has m columns and rank r, iflag=m-r is
// returned. Columns 0:iflag-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. The third argument, E, 
// should be preallocated by the user with m elements, m being
// the number of columns in V. On successful return (iflag>=0),
// e[j] indicates the norm of V(:,j) before the orthogonali-  
// zation step. On exit, V(:,j), j>*iflag has 2-norm 1.
void SUBR(svrr)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, int* iflag);
