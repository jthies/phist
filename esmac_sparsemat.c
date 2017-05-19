#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_constants.h"
#include "esmac_sparsemat.h"

esmac_sparsemat_t * esmac_sparsemat_alloc(int nr, int nc, int ne) {
  esmac_sparsemat_t *sm = malloc(sizeof *sm);
  sm->nr = nr;
  sm->nc = nc;
  sm->ne=0;
  sm->nemax = ne;
  sm->rptr = malloc((nr+1) * sizeof *(sm->rptr));
  sm->cind = malloc( ne    * sizeof *(sm->cind));
  sm->val  = malloc( ne    * sizeof *(sm->val) );
  return sm;
}

void esmac_sparsemat_free(esmac_sparsemat_t *sm) {
  if (sm) {
    if (sm->rptr) {free(sm->rptr);}
    if (sm->cind) {free(sm->cind);}
    if (sm->val ) {free(sm->val );}
  }
}

esmac_sparsemat_t * esmac_sparsemat_from_generator(const esmac_generator_t *gen) {
  
  int info;
  esmac_generator_work_t * my_ws = esmac_generator_alloc(gen, &info);
 
  esmac_idx_t nrow = esmac_generator_query(gen,my_ws,"nrow");
  esmac_idx_t ncol = esmac_generator_query(gen,my_ws,"ncol"); 
  esmac_idx_t mrow = esmac_generator_query(gen,my_ws,"maxnzrow");
  int valtype = (int) esmac_generator_query(gen,my_ws,"valtype");

  esmac_idx_t *cind = malloc(mrow * sizeof *cind);
  double *val = malloc(mrow * sizeof *val);
  
  esmac_idx_t idx;

  // count # non-zeros
  esmac_idx_t ne = 0;
  for (idx=0;idx<nrow;idx++) {
    int k = esmac_generator_row(gen, my_ws, idx, cind, val);
    ne=ne+k;
  }

  esmac_sparsemat_t *sm = esmac_sparsemat_alloc(nrow,ncol,ne);
  sm->valtype = valtype;

  // create matrix
  esmac_idx_t n = 0;
  sm->rptr[0]=0;
  for (idx=0;idx<nrow;idx++) {
    int k = esmac_generator_row(gen, my_ws, idx, cind, val);
    //beware about int <-> esmac_idx_t
    //  memcpy(&(sm->cind[n]), cind, k * sizeof *cind);
    int i;
    for (i=0;i<k;i++) {
      sm->cind[i+n]=cind[i];
    }
    memcpy(&(sm->val [n]), val , k * sizeof *val);
    n=n+k;
    sm->rptr[idx+1]=n;
  }
  if (n != ne) {
    printf("%s: Error.",__func__);
    exit(EXIT_FAILURE);
  }
  sm->ne=ne;

  free(cind);
  free(val);
  esmac_generator_free(my_ws);
  
  return sm;
  
}


int esmac_sparsemat_mvm(const esmac_sparsemat_t *sm, const double *x, double *y, double alpha, double beta, double gamma) {
  int i,j;

  if (gamma != 0.0 && sm->nr != sm->nc) {
    printf("%s: gamma != 0 requires a.nr == a.nc\n!",__func__);
    exit(EXIT_FAILURE);
  }

  // check for early return. TO DO: call BLAS

  if (alpha==0.0 && beta==0.0) {
    for (i=0;i<sm->nr;i++) {y[i]=gamma*x[i];}
    return 0;
  }
  if (alpha==0.0) {
    for (i=0;i<sm->nr;i++) {y[i]=beta*y[i]+gamma*x[i];}
    return 0;
  }
  
  for (i=0;i<sm->nr;i++) {
    //if beta=0 but y[i]=NaN an error might occur because 0*NaN /= 0
    //therefore: check
    if (beta==0.0) {
      y[i]=gamma*x[i];
    } else {
      y[i]=beta*y[i]+gamma*x[i];
    }
    for (j=sm->rptr[i];j<sm->rptr[i+1];j++) {
      y[i]=y[i]+alpha*sm->val[j]*x[sm->cind[j]];
    }
  }
  
  return 0;
}

int esmac_sparsemat_maxrowlength(const esmac_sparsemat_t *sm) {
  int mrwl = 0;
  int i;
  for (i=0;i<sm->nr;i++) {
    mrwl = MAX(mrwl,sm->rptr[i+1]-sm->rptr[i]);
  }
  return mrwl;
}
