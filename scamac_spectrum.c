#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scamac_config.h"
#ifdef SCAMAC_EXTERNAL_CBLAS
#include SCAMAC_EXTERNAL_CBLAS
#else
#include "scamac_cblas.h"
#endif
#ifdef SCAMAC_EXTERNAL_CLAPACK
#include SCAMAC_EXTERNAL_CLAPACK
#else
#include "scamac_clapack.h"
#endif

#include "scamac_spectrum.h"

// Compute entire spectrum of (real or complex, symmetric/hermitian) matrix. Use only for *really* small matrices
ScamacErrorCode scamac_spectrum_real_symmetric(const scamac_sparsemat_st * sm, double **spec) {
  if (!sm || !spec) {
    return SCAMAC_ENULL;
  }

  // obtain matrix

  double *dmat = calloc(sm->nr * sm->nr, sizeof *dmat);  // and is =0
  if (!dmat) {
    return SCAMAC_EMALLOCFAIL;
  }

  ScamacIdx idx;
  for (idx=0; idx<sm->nr; idx++) {
    int i;
    for (i=sm->rptr[idx]; i<sm->rptr[idx+1]; i++) {
      dmat[idx*sm->nr+sm->cind[i]] = sm->val[i];
    }
  }

  // compute spectrum
  char jobz='N';
  char uplo='U';

  int n = sm->nr;

  double *w = malloc(n * sizeof *w);
  if (!w) {
    return SCAMAC_EMALLOCFAIL;
  }
  double *work;
  double wwork;
  int lwork;
  int lapinfo;


  lwork=-1; // workspace query
  dsyev_(&jobz,&uplo,&n,dmat,&n,w,&wwork,&lwork,&lapinfo);
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }
  lwork=wwork+0.1;
  work=malloc(lwork * sizeof *work);
  if (!work) {
    return SCAMAC_EMALLOCFAIL;
  }
  dsyev_(&jobz,&uplo,&n,dmat,&n,w,work,&lwork,&lapinfo);
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }

  free(work);
  free(dmat);

  *spec = w;

  return SCAMAC_EOK;

}

ScamacErrorCode scamac_spectrum_cplx_hermitian(const scamac_sparsemat_st * sm, double **spec) {
  // obtain matrix

  double *dmat = calloc(2* sm->nr * sm->nr, sizeof *dmat);  // and is =0
  if (!dmat) {
    return SCAMAC_EMALLOCFAIL;
  }

  ScamacIdx idx;
  for (idx=0; idx<sm->nr; idx++) {
    int i;
    for (i=sm->rptr[idx]; i<sm->rptr[idx+1]; i++) {
      dmat[2*(idx*sm->nr+sm->cind[i])  ] = sm->val[2*i  ];
      dmat[2*(idx*sm->nr+sm->cind[i])+1] = sm->val[2*i+1];
    }
  }

  // compute spectrum
  char jobz='N';
  char uplo='U';


  int n = sm->nr;

  double *w = malloc(n * sizeof *w);
  if (!w) {
    return SCAMAC_EMALLOCFAIL;
  }
  double complex *work;
  double complex wwork;
  double *rwork = malloc(3*n * sizeof *rwork);
  if (!rwork) {
    return SCAMAC_EMALLOCFAIL;
  }
  int lwork;
  int lapinfo;


  lwork=-1; // workspace query
  zheev_(&jobz,&uplo,&n,(void *) dmat,&n,w,(void *) &wwork,&lwork,rwork,&lapinfo); // nasty VOID instead of complex
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }
  lwork=wwork+0.1;
  work=malloc(lwork * sizeof *work);
  if (!work) {
    return SCAMAC_EMALLOCFAIL;
  }
  zheev_(&jobz,&uplo,&n,(void *) dmat,&n,w,(void *) work,&lwork,rwork,&lapinfo); // nasty VOID instead of complex
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }

  free(work);
  free(dmat);
  free(rwork);

  *spec = w;

  return SCAMAC_EOK;

}

ScamacErrorCode scamac_spectrum_real_general(const scamac_sparsemat_st * sm, double **spec) {
  // obtain matrix

  double *dmat = calloc(sm->nr * sm->nr, sizeof *dmat);  // and is =0
  if (!dmat) {
    return SCAMAC_EMALLOCFAIL;
  }

  ScamacIdx idx;
  for (idx=0; idx<sm->nr; idx++) {
    int i;
    for (i=sm->rptr[idx]; i<sm->rptr[idx+1]; i++) {
      dmat[idx*sm->nr+sm->cind[i]] = sm->val[i];
    }
  }

  // compute spectrum
  char jobvs='N';
  char sort='N';
  int sdim=0;
  int ldvs=1;

  int n = sm->nr;

  double *wr = malloc(n * sizeof *wr);
  if (!wr) {
    return SCAMAC_EMALLOCFAIL;
  }
  double *wi = malloc(n * sizeof *wi);
  if (!wi) {
    return SCAMAC_EMALLOCFAIL;
  }
  double *work;
  double wwork;
  int lwork;
  int lapinfo;

  lwork=-1; // workspace query
  dgees_(&jobvs, &sort, NULL, &n, dmat, &n, &sdim, wr, wi, NULL, &ldvs, &wwork, &lwork, NULL, &lapinfo);
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }
  lwork=wwork+0.1;
  work=malloc(lwork * sizeof *work);
  if (!work) {
    return SCAMAC_EMALLOCFAIL;
  }
  dgees_(&jobvs, &sort, NULL, &n, dmat, &n, &sdim, wr, wi, NULL, &ldvs, work, &lwork, NULL, &lapinfo);
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }

  free(work);
  free(dmat);

  *spec = malloc(2 * n * sizeof **spec);
  if (!spec) {
    return SCAMAC_EMALLOCFAIL;
  }

  cblas_dcopy(n, wr,1, *spec,2);
  cblas_dcopy(n, wi,1, &((*spec)[1]),2);

  free(wr);
  free(wi);

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_spectrum_cplx_general(const scamac_sparsemat_st * sm, double **spec) {
  // obtain matrix

  double *dmat = calloc(2* sm->nr * sm->nr, sizeof *dmat);  // and is =0
  if (!dmat) {
    return SCAMAC_EMALLOCFAIL;
  }

  ScamacIdx idx;
  for (idx=0; idx<sm->nr; idx++) {
    int i;
    for (i=sm->rptr[idx]; i<sm->rptr[idx+1]; i++) {
      dmat[2*(idx*sm->nr+sm->cind[i])  ] = sm->val[2*i  ];
      dmat[2*(idx*sm->nr+sm->cind[i])+1] = sm->val[2*i+1];
    }
  }

  // compute spectrum
  char jobvs='N';
  char sort='N';
  int sdim=0;
  int ldvs=1;

  int n = sm->nr;

  double complex *w = malloc(n * sizeof *w);
  if (!w) {
    return SCAMAC_EMALLOCFAIL;
  }
  double complex *work;
  double complex wwork;
  int lwork;
  double *rwork;
  rwork = malloc(n * sizeof *rwork);
  if (!rwork) {
    return SCAMAC_EMALLOCFAIL;
  }
  int lapinfo;

  lwork=-1; // workspace query
  zgees_(&jobvs, &sort, NULL, &n, (void *) dmat, &n, &sdim, (void *) w, NULL, &ldvs, (void *) &wwork, &lwork, rwork, NULL, &lapinfo);
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }
  lwork=wwork+0.1;
  work=malloc(lwork * sizeof *work);
  if (!work) {
    return SCAMAC_EMALLOCFAIL;
  }
  zgees_(&jobvs, &sort, NULL, &n, (void *) dmat, &n, &sdim, (void *) w, NULL, &ldvs, (void *) work, &lwork, rwork, NULL, &lapinfo);
  if (lapinfo) {
    return (SCAMAC_EFAIL | SCAMAC_EINTERNAL);
  }

  free(work);
  free(dmat);
  free(rwork);

  *spec = (double *) w;

  return SCAMAC_EOK;

}

