// calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
// higher-precision + pivoting + stable low rank approximation
extern "C" void SUBR(prec_cholesky)(_ST_ *__restrict__ a, _ST_ *__restrict__ aC, phist_lidx n, phist_lidx lda, phist_lidx *perm, int *rank, int* iflag)
{
  *iflag=-99;
}

// apply backward substitution with permuted upper triangular matrix to k vectors in col-major storage
extern "C" void SUBR(prec_backwardSubst)(const _ST_ *__restrict__ r, const _ST_ *__restrict__ rC, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        _ST_ *__restrict__ x, _ST_ *__restrict__ xC, phist_lidx k, phist_lidx ldx, int* iflag)
{
  *iflag=-99;
}

// apply forward substitution with permuted transposed upper triangular matrix
extern "C" void SUBR(prec_forwardSubst)(const _ST_ *__restrict__ r, const _ST_ *__restrict__ rC, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        _ST_ *__restrict__ x, _ST_ *__restrict__ xC, phist_lidx k, phist_lidx ldx, int* iflag)
{
  *iflag=-99;
}

// given symmetric A=V'V, compute B s.t. Q=V*B is orthonormal. B is computed in-place,
// if A is found to be numerically rank deficient, the last n-*rank columns of B will be zeroed out
// s.t. Q has exactly rank *rank. In bi, biC we return the inverse of B.
extern "C" void SUBR(prec_qb)(_ST_ *__restrict__ a,  _ST_ *__restrict__ aC, 
                    _ST_ *__restrict__ bi, _ST_ *__restrict__ biC, 
                    phist_lidx n, phist_lidx lda, int *rank, int* iflag)
{
  *iflag=-99;
}

