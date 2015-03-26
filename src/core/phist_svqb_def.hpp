// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182).
// If the input vector has m columns and rank r, iflag=m-r is
// returned. Columns 0:iflag-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. The third argument, E, 
// should be preallocated by the user with m elements, m being
// the number of columns in V. On successful return (iflag>=0),
// e[j] indicates the norm of V(:,j) before the orthogonali-  
// zation step. On exit, V(:,j), j>*iflag has 2-norm 1.
void SUBR(svqb)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, _MT_* D, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
    int m, rank;
    lidx_t ldb;
    _ST_*  B_raw;
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B,&B_raw,&ldb,iflag),*iflag);
    _MT_ Dinv[m]; // inverse sqrt of V'V
    _MT_ E[m], Einv[m]; // sqrt of eigenvalues of (scaled) V'V (and its inverse)
    
    // S=V'V
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),B,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(B,iflag),*iflag);
    // scaling factors: sqrt of inverse diagonal elements
    for (int i=0; i<m; i++)
    {
      _MT_ d=st::real(B_raw[i*ldb+i]);
      // note: diagonal entry must be real
      D[i] = mt::sqrt(d);
      if (mt::abs(d)>mt::eps())
      {
        Dinv[i] = mt::one()/D[i];
      }
      else
      {
        Dinv[i] = mt::zero();
      }
    }
    // scale matrix S^=V'V with ones on diagonal
    for (int i=0; i<m; i++)
    {
      for(int j=0; j<m; j++) 
      {
        B_raw[i*ldb+j] *= Dinv[i]*Dinv[j];
      }
    }

// compute eigenvalues/vectors of scaled B, eigenvalues
// are given in order of ascending magnitude in E, corresponding
// eigenvectors as columns of B
    _CT_ Ec[m];
    TYPE(sdMat_ptr) A = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&A,m,m,NULL,iflag),*iflag);
    lidx_t lda;
    _ST_*  A_raw;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(A,&A_raw,&lda,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(SchurDecomp)(B_raw, ldb, A_raw, lda, m, m, m, LM, mt::eps(), Ec, iflag), *iflag);

PHIST_SOUT(PHIST_INFO,"singular values of W:\n");
for (int i=0;i<m;i++) PHIST_SOUT(PHIST_INFO,"%24.16e\n",sqrt(E[i]));

    // determine rank of input matrix
    MT emax=mt::abs(E[m-1]); 
    rank=m;
    
    if (emax<10*mt::eps())
    {
      rank=0;
    }
    else
    {
      for(int i=0; i<m; i++)
      {
        if (mt::abs(E[i]<10*emax*mt::eps()))
        {
          rank--;
          E[i]=mt::zero();
          Einv[i]=mt::zero();
        }
        else
        {
          Einv[i] = mt::one()/sqrt(E[i]);
          E[i] = sqrt(E[i]);
        }
      }
    }

    // scale the eigenvector matrix, B <- Dinv * B * Einv
    for(int i=0; i<m; i++)
    {
      for(int j=0;j<m;j++)
      {
#ifdef PHIST_SDMATS_ROW_MAJOR
        B_raw[i*ldb+j] *= Dinv[i]*A_raw[i*lda+j]*Einv[j];
#else
        B_raw[j*ldb+i] *= Dinv[i]*A_raw[j*lda+i]*Einv[j];
#endif
      }
    }

    // GPU upload (if this is a GPU process) as we modified B manually
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(B,iflag),*iflag);

    // compute V <- V*B to get an orthogonal V (up to the first (m-rank) columns,
    // which will be exactly zero)
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,B,iflag),*iflag);

// the return value of this function is the rank of the null space of V on entry
*iflag=m-rank;
}
