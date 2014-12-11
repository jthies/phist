// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182).
// If the input vector has m columns and rank r, ierr=m-r is
// returned. Columns 0:ierr-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. It is possible to obtain
// a QR factorization by e.g. computing [q,r]=qr(B), R=r\q, at
// least if V is of full rank.
void SUBR(svqb)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, int* ierr)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    *ierr=0;
    int m, rank;
    lidx_t ldb;
    _ST_*  B_raw;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B,&B_raw,&ldb,ierr),*ierr);
    _MT_ D[m], Dinv[m]; // sqrt of diag of V'V (and its inverse)
    _MT_ E[m], Einv[m]; // sqrt of eigenvalues of (scaled) V'V (and its inverse)
    
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),B,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(B,ierr),*ierr);
    // scaling factors: inverse diagonal elements
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
    // scale matrix V'V
    for(int i=0; i<m; i++)
    {
      for(int j=0; j<m; j++)  
      {
        B_raw[j*ldb+i] *= Dinv[i]*Dinv[j];
      }
    }

// compute eigenvalues/vectors of scaled B, eigenvalues
// are given in order of ascending magnitude in E, corresponding
// eigenvectors as columns of B
#ifdef IS_COMPLEX
    PHIST_CHK_IERR(*ierr=PHIST_LAPACKE(heevd)
        (SDMAT_FLAG, 'V' , 'U', m, (mt::blas_cmplx_t*)B_raw, ldb, E),*ierr);
#else
    PHIST_CHK_IERR(*ierr=PHIST_LAPACKE(syevd)
        (SDMAT_FLAG, 'V' , 'U', m, B_raw, ldb, E),*ierr);
#endif

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
          E[i] = sqrt(E[i]);
          Einv[i] = mt::one()/sqrt(E[i]);
        }
      }
    }

    // scale the eigenvector matrix, B <- Dinv * B * Einv
    for(int i=0; i<m; i++)
    {
      for(int j=0;j<m;j++)
      {
#ifdef PHIST_SDMATS_ROW_MAJOR
        B_raw[i*ldb+j] *= Dinv[i]*Einv[j];
#else
        B_raw[j*ldb+i] *= Dinv[i]*Einv[j];
#endif
      }
    }

    // GPU upload (if this is a GPU process) as we modified B manually
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(B,ierr),*ierr);

    // compute V <- V*B to get an orthogonal V (up to the first (m-rank) columns,
    // which will be exactly zero)
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,B,ierr),*ierr);


// the return value of this function is the rank of the null space of V on entry
*ierr=m-rank;
}
