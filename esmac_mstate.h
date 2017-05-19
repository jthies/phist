#ifndef ESMAC_MSTATE_H
#define ESMAC_MSTATE_H

#include "esmac_constants.h"
#include "esmac_types.h"
#include "esmac_qstate.h"
#include "esmac_dof_fermions.h"
#include "esmac_dof_bosons.h"

typedef struct {
  esmac_idx_t idx;
  double val;
} esmac_idxval_t;

/*
 * multi-state: multiple degrees of freedom
 */

/* the multi-state type */
typedef struct {
  /* number of quantum states in the multi-state */
  int n;
  /* number of states per QS */
  esmac_idx_t *ns;
  /* total number of states */
  esmac_idx_t ns_total;
  /* array of pointers to the respective QS */
  esmac_qstate_t **qs; 
  /* idx of each qs after esmac_mstate_decode */
  esmac_idx_t *idx;

  /* identifier of each dof, if registered (or ESMAC_DOF_NONE) */
  int *which_dof;
  /* and pointer to the dof structure */
  void **dof;

  /* the row */
  int nalloc;
  int nrow;
  esmac_idxval_t *row; 

  /* counter. Needed when adding qstates */
  int *cnt, *maxcnt;
  
} esmac_mstate_t;


/* multi-state routines */

esmac_mstate_t * esmac_mstate_alloc(int n);
void esmac_mstate_free(esmac_mstate_t *ms);

// set QS at position 0 <= pos < n, with ns states per QS
int esmac_mstate_set_qs(esmac_mstate_t *ms, int pos, esmac_idx_t ns, esmac_qstate_t *qs);

// a "high-level" routine: register a DOF directly with the mstate
//int esmac_mstate_register_dof(esmac_mstate_t *ms, int pos, void *dof);
// int esmac_mstate_register_dof_fermions(esmac_mstate_t *ms, int pos, esmac_dof_fermions_t *dof);
int esmac_mstate_register_dof_bosons(esmac_mstate_t *ms, int pos, esmac_dof_bosons_t *dof);


/* multistate routines */

// obtain mstate from idx
// implies esmac_mstate_row_zero
int esmac_mstate_decode(esmac_mstate_t *ms, esmac_idx_t idx);

// individual idx for qs at pos (with 0 <= idx < 
// error for idx < 0
esmac_idx_t esmac_mstate_qs_idx(esmac_mstate_t *ms, int pos);

void esmac_mstate_row_zero(esmac_mstate_t *ms);

// add current active qstates (as manipulated elsewhere) to output row
// if active part of one qs is empty, nothing is added here
// after return, all qstates have become inactive
int esmac_mstate_row_add(esmac_mstate_t *ms, double alpha);

int esmac_mstate_row_add_one_qs(esmac_mstate_t *ms, int pos, double alpha);

/* get one row of a sparse matrix from the multi-state */
int esmac_mstate_row_to_idxval(esmac_mstate_t *ms, int maxrowlength, int keep_zeros, esmac_idx_t *idx, double *val);

#endif /* ESMAC_MSTATE_H */
