#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "esmac_constants.h"
#include "esmac_aux.h"
#include "esmac_qstate.h"


esmac_qstate_t * esmac_qstate_alloc(size_t sz) {

  esmac_qstate_t *qs = malloc(sizeof *qs);

  qs->sourcefst=0;
  qs->nsource=0;
  qs->is_target=0;
  qs->ntarget=0;
  qs->nalloc = 100;
  qs->sz=sz;
  qs->idx=malloc(qs->nalloc * sizeof *(qs->idx));
  qs->val=malloc(qs->nalloc * sizeof *(qs->val));
  qs->rep=malloc(qs->nalloc * qs->sz           );
  qs->wrk=malloc(qs->nalloc * sizeof *(qs->wrk));

  return qs;
}

void esmac_qstate_free(esmac_qstate_t *qs) {
  if (qs) {
    if (qs->idx) {free(qs->idx);}
    if (qs->val) {free(qs->val);}
    if (qs->rep) {free(qs->rep);}
    if (qs->wrk) {free(qs->wrk);}
    free(qs);
  }
}

void esmac_qstate_init(esmac_qstate_t *qs, esmac_idx_t idx, const void * rep) {
  qs->idx[0]=idx;
  qs->val[0]=1.0;
  memcpy(qs->rep, rep, qs->sz);
  esmac_qstate_reset(qs);
}

// to simplify program logic, we copy where, in principle, we need not to
void esmac_qstate_reset(esmac_qstate_t *qs) {
  qs->sourcefst=0;
  qs->nsource=1;
  qs->is_target=0;
  qs->ntarget=0;
}

// copy: let qscp <- qs
//void esmac_qstate_copy(const esmac_qstate_t *qs, esmac_qstate_t *qscp) {
//  printf("%s: not yet implemented\n",__func__);
//  exit(EXIT_FAILURE);
//}


void esmac_qstate_wrap(esmac_qstate_t *qs) {
  if (qs->is_target) {
    esmac_qstate_scrub(qs);
    /*
      int i;
      for (i=0;i<qs->nact;i++) { // a sorry loop
      qs->idx[1+i]=qs->idx[1+qs->nsth+i];
      qs->val[1+i]=qs->val[1+qs->nsth+i];
      memcpy(qs->rep+(1+i)*qs->sz, qs->rep + (1+qs->nsth+i)*qs->sz, qs->sz);
      }
    */
    // if (qs->ntarget && qs->nsource) {
    memmove(&(qs->idx[1]),&(qs->idx[qs->sourcefst+qs->nsource]), qs->ntarget * sizeof *(qs->idx));
    memmove(&(qs->val[1]),&(qs->val[qs->sourcefst+qs->nsource]), qs->ntarget * sizeof *(qs->val));
    memmove(qs->rep + qs->sz, qs->rep + (qs->sourcefst+qs->nsource)*qs->sz, qs->ntarget * qs->sz);
    // }
    qs->sourcefst=1;
    qs->nsource=qs->ntarget;
    qs->is_target=0;
    qs->ntarget=0;
  }
}


void esmac_qstate_scrub(esmac_qstate_t *qs) {
  printf("*** %s: Still to do! ***\n",__func__);
}

void * esmac_qstate_init_get(const esmac_qstate_t *qs, esmac_idx_t *idx) {
  if (idx) {*idx = qs->idx[0];}
  return qs->rep;
}

int esmac_qstate_source_n(const esmac_qstate_t *qs) {
  return qs->nsource;
}

void * esmac_qstate_source_get(const esmac_qstate_t *qs, int k, esmac_idx_t *idx, double *val) {
  if (0 <=k && k<qs->nsource) {
    if (idx) {*idx = qs->idx[k+qs->sourcefst];}
    if (val) {*val = qs->val[k+qs->sourcefst];}
    return qs->rep + (k+qs->sourcefst)*qs->sz;
  } else {
      return NULL;
  }
}

int esmac_qstate_source_get_copy(const esmac_qstate_t *qs, int k, esmac_idx_t *idx, double *val, void *rep) {
  if (0 <=k && k<qs->nsource) {
    if (idx) {*idx = qs->idx[k+qs->sourcefst];}
    if (val) {*val = qs->val[k+qs->sourcefst];}
    if (rep) {memcpy(rep, qs->rep + (k+qs->sourcefst)*qs->sz, qs->sz);}
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}

int esmac_qstate_n(const esmac_qstate_t *qs) {
  if (qs->is_target) {
    return qs->ntarget;
  } else {
    return qs->nsource;
  }
}

void * esmac_qstate_get(const esmac_qstate_t *qs, int k, esmac_idx_t *idx, double *val) {
  if (qs->is_target) {
    if (0 <=k && k<qs->ntarget) {
      if (idx) {*idx = qs->idx[k+qs->nsource+qs->sourcefst];}
      if (val) {*val = qs->val[k+qs->nsource+qs->sourcefst];}
      return qs->rep + (k+qs->nsource+qs->sourcefst)*qs->sz;
    } else {
      return NULL;
    } 
  } else {
    if (0 <=k && k<qs->nsource) {
      if (idx) {*idx = qs->idx[k+qs->sourcefst];}
      if (val) {*val = qs->val[k+qs->sourcefst];}
      return qs->rep + (k+qs->sourcefst)*qs->sz;
    } else {
      return NULL;
    } 
  }
}

int esmac_qstate_target(esmac_qstate_t *qs) {
  if (! qs->is_target) {
    qs->is_target=1;
    qs->ntarget=0;
  }
  return 0;
}

int esmac_qstate_target_zero(esmac_qstate_t *qs) {
  qs->is_target=1;
  qs->ntarget=0;
  return 0;
}

// add individual basis state to "active" qstate 
int esmac_qstate_target_add(esmac_qstate_t *qs, esmac_idx_t idx, double val, const void *rep) {
  qs->is_target=1;
  if (val != 0.0) {
    if (esmac_increase_n_somewhat(qs->nsource+qs->ntarget+2)>qs->nalloc) {
      qs->nalloc = esmac_increase_n_somewhat(qs->nsource+qs->ntarget+2);
      qs->idx = realloc(qs->idx, qs->nalloc * sizeof *(qs->idx) );
      qs->val = realloc(qs->val, qs->nalloc * sizeof *(qs->val) );
      qs->rep = realloc(qs->rep, qs->nalloc * qs->sz            );
    }
    int j = qs->sourcefst+qs->nsource+qs->ntarget;
    qs->idx[j] = idx;
    qs->val[j] = val;
    memcpy(qs->rep + j*qs->sz, rep, qs->sz);
    (qs->ntarget)++;
  }
  return 0;
}

// multiply by scalar: qs -> alpha * qs 
void esmac_qstate_target_scale(esmac_qstate_t *qs, double alpha) {
  qs->is_target=1;
  if (alpha==0.0) {
    qs->ntarget=0;
  } else if (alpha != 1.0) {
    int i;
    for (i=0;i<qs->ntarget;i++) {
      int j = qs->sourcefst+qs->nsource+i;
      qs->val[j] = alpha * qs->val[j];
    }
  }
}


int esmac_qstate_op_identity(esmac_qstate_t *qs, double alpha) {
  qs->is_target=1;
  qs->ntarget=0;
  if (alpha != 0.0) {
    qs->ntarget=qs->nsource;
    if (esmac_increase_n_somewhat(qs->nsource+qs->ntarget+1)>qs->nalloc) {
      qs->nalloc = esmac_increase_n_somewhat(qs->nsource+qs->ntarget+1);
      qs->idx = realloc(qs->idx, qs->nalloc * sizeof *(qs->idx) );
      qs->val = realloc(qs->val, qs->nalloc * sizeof *(qs->val) );
      qs->rep = realloc(qs->rep, qs->nalloc * qs->sz            );
    }
    memcpy(&(qs->idx[qs->sourcefst+qs->nsource]),&(qs->idx[qs->sourcefst]),qs->ntarget * sizeof *(qs->idx));
    memcpy(&(qs->val[qs->sourcefst+qs->nsource]),&(qs->val[qs->sourcefst]),qs->ntarget * sizeof *(qs->val));
    memcpy(qs->rep + (qs->sourcefst+qs->nsource)*qs->sz, qs->rep + qs->sourcefst*qs->sz,  qs->ntarget * qs->sz);
    if (alpha != 1.0) {
      int i;
      for (i=qs->sourcefst+qs->nsource;i<qs->sourcefst+qs->nsource+qs->ntarget;i++) {
	qs->val[i] = alpha * qs->val[i];
      }
    }
  }
  return 0;
}
