#ifndef ESMAC_DOF_FERMIONS_H
#define ESMAC_DOF_FERMIONS_H

#include "esmac_types.h"

/* quantum mechanical degree of freedom: fermions */
 
typedef struct {
  /* number of sites */
  int n_sites;
  /* number of fermions */
  int n_fermions;

  /* number of states */
  esmac_idx_t ns;
  
  /* count information */
  esmac_idx_t *cnt;
 
} esmac_dof_fermions_t;

//typedef unsigned char * esmac_rep_fermions_t; // "bit" would suffice
typedef int esmac_rep_fermions_t;

esmac_dof_fermions_t * esmac_dof_fermions_alloc(int n_sites, int n_fermions);
void esmac_dof_fermions_free(esmac_dof_fermions_t * dof);

esmac_idx_t esmac_dof_fermions_ns(const esmac_dof_fermions_t * dof);

esmac_rep_fermions_t * esmac_rep_fermions_alloc(const esmac_dof_fermions_t *dof);
void esmac_rep_fermions_free(esmac_rep_fermions_t * rep);

void esmac_rep_fermions_copy(const esmac_dof_fermions_t *dof, const esmac_rep_fermions_t * rep, esmac_rep_fermions_t * repcpy);

int esmac_fermions_decode(const esmac_dof_fermions_t * dof, esmac_idx_t idx, esmac_rep_fermions_t *rep);
esmac_idx_t esmac_fermions_encode(const esmac_dof_fermions_t * dof, const esmac_rep_fermions_t *rep);

/* operators */
double esmac_op_fermions_cdc  (esmac_dof_fermions_t * dof, const esmac_rep_fermions_t *rep, int ic);
double esmac_op_fermions_nn   (esmac_dof_fermions_t * dof, const esmac_rep_fermions_t *rep, int i, int j);

double esmac_op_fermions_cdicj(esmac_dof_fermions_t * dof, esmac_rep_fermions_t *rep, int i, int j);
double esmac_op_fermions_hop  (esmac_dof_fermions_t * dof, esmac_rep_fermions_t *rep, int i, int j);


#endif /* ESMAC_DOF_FERMIONS_H */
