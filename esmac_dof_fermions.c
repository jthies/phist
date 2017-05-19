#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_constants.h"
#include "esmac_lut.h"
#include "esmac_dof_fermions.h"

esmac_dof_fermions_t * esmac_dof_fermions_alloc(int n_sites, int n_fermions) {
  esmac_dof_fermions_t * mydof = malloc(sizeof *mydof);

  mydof->n_sites=n_sites;
  mydof->n_fermions=n_fermions;

  //mydof->ns = esmac_lut_construct(0, n_sites, 1, n_fermions, &(mydof->cnt));
  mydof->ns = esmac_lut_onezero_construct(n_sites, n_fermions, &(mydof->cnt));
  
  return mydof;
}

void esmac_dof_fermions_free(esmac_dof_fermions_t * dof) {
  if (dof) {
    if (dof->cnt) {free(dof->cnt);}
    free(dof);
  }
}

esmac_idx_t esmac_dof_fermions_ns(const esmac_dof_fermions_t * dof) {
  return dof->ns;
}

esmac_rep_fermions_t * esmac_rep_fermions_alloc(const esmac_dof_fermions_t *dof) {
  esmac_rep_fermions_t * myrep = malloc(dof->n_sites * sizeof *myrep);
  return myrep;
}

void esmac_rep_fermions_free(esmac_rep_fermions_t * rep) {
  if (rep) {free(rep);}
}

void esmac_rep_fermions_copy(const esmac_dof_fermions_t *dof, const esmac_rep_fermions_t * rep, esmac_rep_fermions_t * repcpy) {
  if (rep && repcpy) {
    memcpy(repcpy, rep, dof->n_sites * sizeof *rep );
  }
}

int esmac_fermions_decode(const esmac_dof_fermions_t * dof, esmac_idx_t idx, esmac_rep_fermions_t *rep) {
  if (0 <= idx && idx < dof->ns) {
   // esmac_lut_decode(0, dof->n_sites, dof->n_fermions, dof->cnt, idx, dof->rep);
    esmac_lut_onezero_decode(dof->n_sites, dof->n_fermions, dof->cnt, idx, rep);
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}

esmac_idx_t esmac_fermions_encode(const esmac_dof_fermions_t * dof, const esmac_rep_fermions_t *rep) {
  return esmac_lut_onezero_encode(dof->n_sites, dof->n_fermions, dof->cnt, rep);
}

/*
static esmac_idx_t esmac_dof_fermions_encode(const esmac_dof_fermions_t * dof, const void* rep) {
  //return esmac_lut_encode(0, dof->n_sites, dof->n_fermions, dof->cnt, rep);
  return esmac_lut_onezero_encode(dof->n_sites, dof->n_fermions, dof->cnt, rep);
}
*/

/* manipulation of representations */

// c^+_i c_i
static int fermions_cdc(int i, const int *x) {
  return x[i];
}

// x -> scalar * c^+_i c_j x
static int fermions_cdicj(int i, int j, int *x) {
  if (x[i]==1 || x[j]==0) {
    return 0;
  } else {
    // fermionic sign
    int s1,s2;
    s1=0; s2=0;
    int k;
    for (k=0;k<i;k++) {s1=s1+x[k];}
    x[j]=0;
    for (k=0;k<j;k++) {s2=s2+x[k];}
    x[i]=1;
    
    if ((s1*s2)%2 == 0) {
      return 1;
    } else {
      return -1;
    }
  }
}

// n_i n_j
static int fermions_nn(int i, int j, const int *x) {
  if (x[i] && x[j]) {
    return 1;
  } else {
    return 0;
  }
}

// x -> scalar * (c^+_i c_j + c^+_j c_i) x
static int fermions_hop(int i, int j, int *x) {
  if (i==j) {
    if (x[i]) {
      return 2;
    } else {
      return 0;
    }
  } else if ( (x[i]==0 && x[j]==0) || (x[i]==1 && x[j]==1) ) {
    return 0;
  } else {	
    int il,ir;
    if (i<j) {
      il=i;
      ir=j;
    } else {
      il=j;
      ir=i;
    }
    int fsgn,k;
    fsgn=0;
    for (k=il+1;k<ir;k++) {
      if (x[k]) fsgn=1-fsgn;
    }
    x[i]=1-x[i];
    x[j]=1-x[j];
    if (fsgn) {
      return -1;
    } else {
      return 1;
    } 
  }
}


double esmac_op_fermions_cdc  (esmac_dof_fermions_t * dof, const esmac_rep_fermions_t *rep, int ic) {
  if (0 <= ic && ic < dof->n_sites) {
    return (double) fermions_cdc(ic,rep);
  } else {
    return 0.0;
  }
}

double esmac_op_fermions_nn   (esmac_dof_fermions_t * dof, const esmac_rep_fermions_t *rep, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    return (double) fermions_nn(i,j,rep);
  } else {
    return 0.0;
  }
}

double esmac_op_fermions_cdicj  (esmac_dof_fermions_t * dof, esmac_rep_fermions_t *rep, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    return (double) fermions_cdicj(i,j,rep);
  } else {
    return 0.0;
  }
}

double esmac_op_fermions_hop  (esmac_dof_fermions_t * dof, esmac_rep_fermions_t *rep, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    return (double) fermions_hop(i,j,rep);
  } else {
    return 0.0;
  }
}
