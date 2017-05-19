#ifndef ESMAC_DOF_BOSONS_H
#define ESMAC_DOF_BOSONS_H

#include "esmac_qstate.h"

/* quantum mechanical degree of freedom: bosons */
 
typedef struct {
  /* number of sites */
  int n_sites;
  /* number of bosons */
  int n_bosons;
  /* maximal number of bosons per site */
  int n_bosons_per_site;

  /* number of states */
  esmac_idx_t ns;
  /* size of representation vector */
  size_t repsize;
  /* representation to work with */
  int *rep;
  
  /* count information */
  esmac_idx_t *cnt;

  /* the associated quantum state */
  esmac_qstate_t *qs;
 
} esmac_dof_bosons_t;


esmac_dof_bosons_t * esmac_dof_bosons_alloc(int n_sites, int n_bosons, int n_bosons_per_site);
void esmac_dof_bosons_free(esmac_dof_bosons_t * dof);

esmac_idx_t esmac_dof_bosons_ns(const esmac_dof_bosons_t * dof);
esmac_qstate_t * esmac_dof_bosons_qs(const esmac_dof_bosons_t * dof);

int esmac_dof_bosons_init(esmac_dof_bosons_t * dof, esmac_idx_t idx);

/* operators */
int esmac_op_bosons_bdb  (esmac_dof_bosons_t * dof, double alpha, int ic);
int esmac_op_bosons_bdibj(esmac_dof_bosons_t * dof, double alpha, int i, int j);
int esmac_op_bosons_nn   (esmac_dof_bosons_t * dof, double alpha, int i, int j);



#endif /* ESMAC_DOF_BOSONS_H */
