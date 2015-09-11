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

//! symmetric eigenvalue decomposition in simulated quad precision (using function
//! Rsyev, cf. lapack routine dsyev)
extern "C" void phist_Drsyev(int n, double *restrict a, double *restrict aC, int lda,
                       double *restrict w, double *restrict wC, int *iflag)
{
#if defined(TESTING) && (PHIST_OUTLEV>=PHIST_TRACE)
  printf("Entering %s\n", __FUNCTION__);
#endif

    mpackint lwork, liwork;

    dd_real A[n*n];
    dd_real W[n];

//copy A matrix to GMP data structure
for (int i=0; i<n; i++)
{
  for (int j=0; j<n; j++)
  {
    A[i+j*n].x[0] = a[i+j*lda];
    A[i+j*n].x[1] = aC[i+j*lda];
  }
}
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

    for (int i=0; i<n; i++)
    {
      w[i]=W[i].x[0];
      wC[i]=W[i].x[1];
    }

//copy A matrix to GMP data structure
for (int i=0; i<n; i++)
{
  for (int j=0; j<n; j++)
  {
    a [i+j*lda]=A[i+j*n].x[0];
    aC[i+j*lda]=A[i+j*n].x[1];
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


/* PHIST_HAVE_MPACK_QD */
#endif
