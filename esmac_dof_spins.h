#ifndef ESMAC_DOF_SPINS_H
#define ESMAC_DOF_SPINS_H

#include "esmac_types.h"

/* quantum mechanical degree of freedom: spins 1/2 */
/* either with J_z = const., or without J_z constraint */

typedef struct {
  /* number of sites */
  int n_sites;
  /* number of "up" spins 
     negative value indicates "no constraint" */
  int n_up;

  /* number of states */
  esmac_idx_t ns;
  
  /* count information */
  esmac_idx_t *cnt;
 
} esmac_dof_spins_t;



//typedef unsigned char * esmac_rep_fermions_t; // "bit" would suffice
typedef int esmac_rep_spins_t;

esmac_dof_spins_t * esmac_dof_spins_alloc(int n_sites);
esmac_dof_spins_t * esmac_dof_spinssz_alloc(int n_sites, int n_up);
void esmac_dof_spins_free(esmac_dof_spins_t * dof);


esmac_idx_t esmac_dof_spins_ns(const esmac_dof_spins_t * dof);

esmac_rep_spins_t * esmac_rep_spins_alloc(const esmac_dof_spins_t *dof);
void esmac_rep_spins_free(esmac_rep_spins_t * rep);

void esmac_rep_spins_copy(const esmac_dof_spins_t *dof, const esmac_rep_spins_t * rep, esmac_rep_spins_t * repcpy);

int esmac_spins_decode(const esmac_dof_spins_t * dof, esmac_idx_t idx, esmac_rep_spins_t *rep);
esmac_idx_t esmac_spins_encode(const esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep);


/* operators */
double esmac_op_spins_Sz   (esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep, int i);
double esmac_op_spins_SzSz (esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep, int i, int j);
double esmac_op_spins_SzSum(esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep);


double esmac_op_spins_Sp   (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i);
double esmac_op_spins_Sm   (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i);
double esmac_op_spins_Sx   (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i);

double esmac_op_spins_SpSp (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j);
double esmac_op_spins_SpSm (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j);
double esmac_op_spins_SmSm (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j);
double esmac_op_spins_SxSx (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j);
/* Sx Sx + Sy Sy = (S+ S- + S- S+)/2 */
double esmac_op_spins_Sflip(esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j);

#endif /* ESMAC_DOF_SPINS_H */
