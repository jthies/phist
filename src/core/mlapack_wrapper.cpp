/*! \file sdmat_syev_prec.cpp
 * interface to multiple-precision lapack routine Rsyev
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"
#include "prec_helpers.h"
#include "phist_defs.h"

#ifdef PHIST_HAVE_MPACK_QD

#include <cstdio>

#include "qd/dd_real.h"
#include "mpack/dd_complex.h"
#include "mpack/mpack_config.h"
#include "mpack/mlapack_dd.h"

#if( PHIST_OUTLEV >= PHIST_DEBUG )
void printmat(int N, int M, dd_real * A, int LDA)
{
    dd_real mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            mtmp = A[i + j * LDA];
            printf("  %24.16e", mtmp.x[0]+mtmp.x[1]);
            if (j < M - 1)
                printf(", ");
        }
        if (i < N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}
#endif

// copy sdMat in (a,aC) to GMP data structure
#define PHIST_TO_QD(_a,_aC,_lda,_nrows,_ncols,_A) \
{\
for (int i=0; i<(_ncols); i++)\
{\
  for (int j=0; j<(_nrows); j++)\
  {\
    (_A)[i+j*(_nrows)].x[0] = (_a) [i+j*(_lda)];\
    (_A)[i+j*(_ncols)].x[1] = (_aC)[i+j*(_lda)];\
  }\
}\
}

// copy sdMat in (a,aC) to GMP data structure
#define QD_TO_PHIST(_a,_aC,_lda,_nrows,_ncols,_A) \
{\
for (int i=0; i<(_ncols); i++)\
{\
  for (int j=0; j<(_nrows); j++)\
  {\
    (_a) [i+j*(_lda)] = (_A)[i+j*(_nrows)].x[0]; \
    (_aC)[i+j*(_lda)] = (_A)[i+j*(_nrows)].x[1]; \ 
  }\
}\
}

//! symmetric eigenvalue decomposition in simulated quad precision (using function
//! Rsyev, cf. lapack routine dsyev)
//! in contrast to dsyev we return the eigenvalues in reversed order!
extern "C" void phist_Drsyev(int n, double *restrict a, double *restrict aC, int lda,
                       double *restrict w, double *restrict wC, int *iflag)
{
#if defined(PHIST_TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

    mpackint lwork, liwork;

    dd_real A[n*n];
    dd_real W[n];
    
    PHIST_TO_QD(a,aC,lda,n,n,A);

#if( PHIST_OUTLEV >= PHIST_DEBUG )
    printf("%%in Rsyev,\n A =");
    printmat(n, n, A, n);
    printf("\n");
#endif
//work space query
    lwork = (mpackint)(-1);
    mpackint info;
    dd_real tmp_work;
    mpackint nn=(mpackint)n;
    Rsyev("V", "U", nn, A, nn, W, &tmp_work, lwork, &info);
    lwork = std::max((mpackint)tmp_work.x[0],(mpackint)1);
      
    dd_real *work=new dd_real[lwork];

//compute eigenvalues
    Rsyev("V", "U", nn, A, nn, W, work, lwork, &info);

    delete [] work;

// set eigenvalues in reversed order
    for (int i=0; i<n; i++)
    {
      w[i]=W[n-i-1].x[0];
      wC[i]=W[n-i-1].x[1];
    }

//copy A matrix from GMP data structure in reversed order
for (int i=0; i<n; i++)
{
  for (int j=0; j<n; j++)
  {
    a [i+j*lda]=A[i+(n-j-1)*n].x[0];
    aC[i+j*lda]=A[i+(n-j-1)*n].x[1];
  }
}

#if( PHIST_OUTLEV >= PHIST_DEBUG )
    printf("%%in Rsyev, eigenvalues \n");
    printf("w =");
    printmat(n, 1, W, 1);
    printf("\n");
    printf("%%in Rsyev, eigenvecs \n");
    printf("U =");
    printmat(n, n, A, n);
    printf("\n");
#endif
}

void phist_Drgesvd(const char *jobu, const char *jobvt, int m, int n,
            double *restrict a, double *restrict aC, int lda, double *restrict s, double *restrict sC,
            double *restrict u, double *restrict uC, int ldu, double *restrict vt, double *restrict vtC, 
            int ldvt, int *iflag)
{
  dd_real A[n*m], S[std::min(n,m)],U[m*m],Vt[n*n];  
  PHIST_TO_QD(a,aC,lda,m,n,A);
  // create work array
  dd_real* work=NULL;
  mpackint lwork=-1;
  dd_real tmp_work;
  Rgesvd(jobu,jobvt,m,n,A,m,S,U,m,Vt,n,&tmp_work,lwork,iflag);
  lwork=(mpackint)tmp_work.x[0];
  Rgesvd(jobu,jobvt,m,n,A,m,S,U,m,Vt,n,iflag);
  delete [] work;
  QD_TO_PHIST( s, sC, m, m, 1, S);
  QD_TO_PHIST( u, uC, m, m, m, U);
  QD_TO_PHIST(vt,vtC, n, n, n, Vt);
}


/* PHIST_HAVE_MPACK_QD */
#endif
