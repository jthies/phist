// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182).
// If the input vector has m columns and rank r, ierr=m-r is
// returned. Columns 0:ierr-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. It is possible to obtain
// a QR factorization by e.g. computing [q,r]=qr(B), R=r\q, at
// least if V is of full rank.
void SUBR(svqb)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, int* ierr)
{
    ENTER_FCN(__FUNCTION__);
    *ierr=0;
    int m, rank;
    lidx_t ldb;
    _ST_*  b;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B,&b,&ldb,ierr),*ierr);
    _MT_ D[m], Dinv[m]; // sqrt of diag of V'V (and its inverse)
    _MT_ E[m], Einv[m]; // sqrt of eigenvalues of (scaled) V'V (and its inverse)
    
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),B,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(B,ierr),*ierr);
    // scaling factors: inverse diagonal elements
    for (int i=0; i<n; i++) 
    {
      _MT_ d=(_MT_)b[i*ldb+i];
      // note: diagonal entry must be real
      D[i] = mt::sqrt(d);
      if (mt::abs(d)>mt::eps())
      {
        Dinv[i] = mt::one()/mt::sqrt(b[i*ldb+i]);
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
        b[j*ldr+i] *= Dinv[i]*Dinv[j];
      }
    }
// compute eigenvalues/vectors of scaled B, eigenvalues
// are given in order of ascending magnitude in E, corresponding
// eigenvectors as columns of B
#ifdef IS_COMPLEX
    PHIST_CHK_IERR(*ierr=LAPACKE_PREFIX(heevd)
        (SDMAT_FLAG, 'V' , 'U', n, b, ldb, E),*ierr);
#else
    PHIST_CHK_IERR(*ierr=LAPACKE_PREFIX(syevd)
        (SDMAT_FLAG, 'V' , 'U', n, b, ldb, E),*ierr);
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
        if (mt::abs(eigs[i]<10*emax*mt::eps()))
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
    for( j=0;j<n;j++)  
    {
#ifdef PHIST_SDMATS_ROW_MAJOR
      b[i*ldb+j] *= Dinv[i]*E[j];
#else
      b[j*ldb+i] *= Dinv[i]*E[j];
#endif
      }
    }
    // compute V <- V*B to get an orthogonal V (up to the first (m-rank) columns,
    // which will be exactly zero)
    PHIST_ICHK_IERR(SUBR(sdMat_to_device)(B,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_times_sdmat_in_place)(st::one(),V,B,ierr),*ierr);


// the return value of this function is the rank of the null space of V on entry
*ierr=m-rank;
}
