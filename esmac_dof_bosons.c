#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "esmac_constants.h"
#include "esmac_lut.h"
#include "esmac_dof_bosons.h"

esmac_dof_bosons_t * esmac_dof_bosons_alloc(int n_sites, int n_bosons, int n_bosons_per_site) {
  esmac_dof_bosons_t * mydof = malloc(sizeof *mydof);

  mydof->n_sites=n_sites;
  mydof->n_bosons=n_bosons;
  if (n_bosons_per_site < 0) {n_bosons_per_site = n_bosons;}
  mydof->n_bosons_per_site=n_bosons_per_site;

  mydof->ns = esmac_lut_construct(0, n_sites, n_bosons_per_site, n_bosons, &(mydof->cnt));
    
  mydof->repsize = n_sites * sizeof(int);
  mydof->rep = malloc(n_sites * sizeof *(mydof->rep));

  mydof->qs = esmac_qstate_alloc(mydof->repsize);
  
  return mydof;
}
  
void esmac_dof_bosons_free(esmac_dof_bosons_t * dof) {
  if (dof) {
    if (dof->rep) {free(dof->rep);}
    if (dof->cnt) {free(dof->cnt);}
    if (dof->qs) {esmac_qstate_free(dof->qs);}
    free(dof);
  }
}

esmac_idx_t esmac_dof_bosons_ns(const esmac_dof_bosons_t * dof) {
  return dof->ns;
}

esmac_qstate_t * esmac_dof_bosons_qs(const esmac_dof_bosons_t * dof) {
  return dof->qs;
}

int esmac_dof_bosons_init(esmac_dof_bosons_t * dof, esmac_idx_t idx) {
  if (0 <= idx && idx < dof->ns) {
    esmac_lut_decode(0, dof->n_sites, dof->n_bosons, dof->cnt, idx, dof->rep);
    esmac_qstate_init(dof->qs, idx, dof->rep);
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}


static esmac_idx_t esmac_dof_bosons_encode(const esmac_dof_bosons_t * dof, const void* rep) {
  int * my_rep = (int *) rep;
  int i;
  for (i=0;i<dof->n_sites;i++) {
    if (my_rep[i] > dof->n_bosons_per_site) {
      return -1; // states does not exist
    }
  }
  return esmac_lut_encode(0, dof->n_sites, dof->n_bosons, dof->cnt, rep);
}


/* manipulation of representations */

// b^+_i b_i
static int bosons_bdb(int i, const int *x) {
  return x[i];
}

// x -> scalar * b^+_i b_j x
static double bosons_bdibj(int i, int j, int *x) {
  if (x[j]==0) {
    return 0;
  } else {
    double val = sqrt((double) (x[i]+1)*x[j]);
    x[i]++;
    x[j]--;
    return val;
  }
}

// n_i n_j
static int bosons_nn(int i, int j, const int *x) {
  return x[i]*x[j];
}


int esmac_op_bosons_bdb  (esmac_dof_bosons_t * dof, double alpha, int ic) {
  if (0 <= ic && ic < dof->n_sites) {
    esmac_qstate_target(dof->qs);
    if (alpha!=0.0) {
      int k;
      for (k=0;k<esmac_qstate_source_n(dof->qs);k++) {
        esmac_idx_t idx;
        double val;
        int *rep;
        rep = esmac_qstate_source_get(dof->qs, k, &idx, &val);
        int ires = bosons_bdb(ic,rep);
        if (ires) { // non-zero result == 1
          esmac_qstate_target_add(dof->qs, idx, alpha*( (double) ires)*val, rep);
        }
      }
    }
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}

int esmac_op_bosons_bdibj(esmac_dof_bosons_t * dof, double alpha, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    esmac_qstate_target(dof->qs);
    if (alpha!=0.0) {
      int k;
      for (k=0;k<esmac_qstate_source_n(dof->qs);k++) {
        esmac_idx_t idx;
        double val;
        esmac_qstate_source_get_copy(dof->qs, k, &idx, &val, dof->rep);
        double res = bosons_bdibj(i,j,dof->rep);
        if (res != 0.0) { // non-zero result
          esmac_idx_t new_idx = esmac_dof_bosons_encode(dof, dof->rep);
          if (idx>=0) {
            esmac_qstate_target_add(dof->qs, new_idx, alpha * res * val, dof->rep);
          }
        }
      }
    }
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}


int esmac_op_bosons_nn   (esmac_dof_bosons_t * dof, double alpha, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    esmac_qstate_target(dof->qs);
    if (alpha!=0.0) {
      int k;
      for (k=0;k<esmac_qstate_source_n(dof->qs);k++) {
        esmac_idx_t idx;
        double val;
        int *rep;
        rep = esmac_qstate_source_get(dof->qs, k, &idx, &val);
        int ires = bosons_nn(i,j,rep);
        if (ires) { // non-zero result == 1
          esmac_qstate_target_add(dof->qs, idx, alpha* ((double) ires) * val, rep);
        }
      }
    }
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}



