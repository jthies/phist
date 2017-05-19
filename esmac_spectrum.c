#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_spectrum.h"

extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

// Compute entire spectrum of (real, symmetric) matrix. Use only for *really* small matrices
int esmac_spectrum(const esmac_generator_t * gen, double **spec) {
  // obtain matrix

  int info;
  esmac_generator_work_t * my_ws = esmac_generator_alloc(gen, &info);
 
  esmac_idx_t nrow = esmac_generator_query(gen,my_ws,"nrow");
  esmac_idx_t ncol = esmac_generator_query(gen,my_ws,"ncol"); 
  esmac_idx_t mrow = esmac_generator_query(gen,my_ws,"maxnzrow");
  if (ncol != nrow) {
    printf("%s: Requires square matrix. Abort.\n",__func__);
    exit(EXIT_FAILURE);
  }

  esmac_idx_t *cind = malloc(mrow * sizeof *cind);
  double *val = malloc(mrow * sizeof *val);
  
  double *dmat = calloc(nrow * nrow , sizeof *dmat); // and is =0
    
  esmac_idx_t idx;
  for (idx=0;idx<nrow;idx++) {
    int n = esmac_generator_row(gen, my_ws, idx, cind, val);
  
    int i;
    for (i=0;i<n;i++) {
      dmat[idx*nrow+cind[i]] = val[i];
    }
  }
  free(cind);
  free(val);
  esmac_generator_free(my_ws);
  
  
  // compute spectrum
  char jobz='N';
  char uplo='U';

  double *w = malloc(nrow * sizeof *w);
  double *work;
  double wwork;
  int lwork;
  int lapinfo;
  
  int n = nrow;


  lwork=-1; // workspace query
  dsyev_(&jobz,&uplo,&n,dmat,&n,w,&wwork,&lwork,&lapinfo);
  lwork=wwork+0.1;
  work=malloc(lwork * sizeof *work);
  dsyev_(&jobz,&uplo,&n,dmat,&n,w,work,&lwork,&lapinfo);

  free(work);
  free(dmat);
  
  *spec = w;
  
  return n;

}

  
  // Compute entire spectrum of (real, symmetric) matrix. Use only for *really* small matrices
int esmac_spectrum_mat(const esmac_sparsemat_t * sm, double **spec) {
  // obtain matrix

  double *dmat = calloc(sm->nr * sm->nr , sizeof *dmat); // and is =0
    
  esmac_idx_t idx;
  for (idx=0;idx<sm->nr;idx++) {
    int i;
    for (i=sm->rptr[idx];i<sm->rptr[idx+1];i++) {
      dmat[idx*sm->nr+sm->cind[i]] = sm->val[i];
    }
  }
  
  // compute spectrum
  char jobz='N';
  char uplo='U';


  int n = sm->nr;

  double *w = malloc(n * sizeof *w);
  double *work;
  double wwork;
  int lwork;
  int lapinfo;
  

  lwork=-1; // workspace query
  dsyev_(&jobz,&uplo,&n,dmat,&n,w,&wwork,&lwork,&lapinfo);
  lwork=wwork+0.1;
  work=malloc(lwork * sizeof *work);
  dsyev_(&jobz,&uplo,&n,dmat,&n,w,work,&lwork,&lapinfo);

  free(work);
  free(dmat);
  
  *spec = w;
  
  return n;

}

