#ifndef ESMAC_QSTATE_H
#define ESMAC_QSTATE_H

#include <stddef.h> // size_t

#include "esmac_types.h"

/* a state of a single quantum mechanical degree of freedom */
typedef struct {
  // first state of current source state (is either 0 or 1)
  int sourcefst;
  // number of elements in source state
  int nsource;
  // has an operator been applied to the state?
  int is_target;
  // number of elements in target state
  int ntarget;
  // number of allocated elements
  int nalloc;
  // size of representations
  size_t sz;
  // the indices
  esmac_idx_t * idx;
  // the amplitudes
  double * val;
  // the representations, each of size sz
  char * rep;
  // work array for "scrub"
  esmac_idx_t *wrk;
} esmac_qstate_t;

esmac_qstate_t * esmac_qstate_alloc(size_t sz);
void esmac_qstate_free(esmac_qstate_t *qs);

// initialize with one basis state (idx, rep)
void esmac_qstate_init(esmac_qstate_t *qs, esmac_idx_t idx, const void * rep);
// reset to initialized state
void esmac_qstate_reset(esmac_qstate_t *qs);

// copy: let qscp <- qs
//void esmac_qstate_copy(const esmac_qstate_t *qs, esmac_qstate_t *qscp);


// add two states: qs -> qs + alpha qsx
//void esmac_qstate_add_two_qs(esmac_qstate_t *qs, double alpha, esmac_qstate_t *qsx);

// make target state to source state
void esmac_qstate_wrap(esmac_qstate_t *qs);

// remove doublets etc. from active part. Is implied by esmac_qstate_stash, before stashing
void esmac_qstate_scrub(esmac_qstate_t *qs);

// get individual INIT basis state
void * esmac_qstate_init_get(const esmac_qstate_t *qs, esmac_idx_t *idx);

// number of individual basis states in "stashed" qstate
int esmac_qstate_source_n(const esmac_qstate_t *qs);
// return individual basis state, 0<=k < esmac_qstate_length
void * esmac_qstate_source_get(const esmac_qstate_t *qs, int k, esmac_idx_t *idx, double *val);
// get copy of individual basis state
int esmac_qstate_source_get_copy(const esmac_qstate_t *qs, int k, esmac_idx_t *idx, double *val, void *rep);

// number of individual basis states in source or target qstate (depending on whether state is targeted)
int esmac_qstate_n(const esmac_qstate_t *qs);
// return individual basis state, 0<=k < esmac_qstate_length
void * esmac_qstate_get(const esmac_qstate_t *qs, int k, esmac_idx_t *idx, double *val);


// "activate" qstate (implies esmac_qstate_active_zero, if it hasn't been a target before)
int esmac_qstate_target(esmac_qstate_t *qs);
// set active state to zero (implies esmac_qstate_activate)
int esmac_qstate_target_zero(esmac_qstate_t *qs);
// add individual basis state to "active" qstate 
int esmac_qstate_target_add(esmac_qstate_t *qs, esmac_idx_t idx, double val, const void *rep);
// multiply active state by scalar: qs -> alpha * qs 
void esmac_qstate_target_scale(esmac_qstate_t *qs, double alpha);

// basically, target = alpha * source
int esmac_qstate_op_identity(esmac_qstate_t *qs, double alpha);

#endif /* ESMAC_QSTATE_H */
