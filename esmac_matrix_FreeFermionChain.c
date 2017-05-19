#include <stdlib.h>

#include "esmac_constants.h"
#include "esmac_dof_fermions.h"
#include "esmac_mstate.h"
#include "esmac_matrix_FreeFermionChain.h"


esmac_matrix_FreeFermionChain_work_t * esmac_matrix_FreeFermionChain_alloc(const esmac_matrix_FreeFermionChain_params_t * par, int * info) {
  esmac_matrix_FreeFermionChain_work_t * ws = malloc(sizeof *ws);
  if (par->n_species != 0) {
    ws->ndof = par->n_species;
    ws->midx = esmac_multidx_alloc(par->n_species);
    ws->row = esmac_sparserow_real_alloc();
    ws->dof = malloc(par->n_species * sizeof *(ws->dof));
    ws->rep = malloc(par->n_species * sizeof *(ws->rep));
    ws->repinit = malloc(par->n_species * sizeof *(ws->repinit));
    ws->ns = 1;
    int k;
    for (k=0;k<par->n_species;k++) {
      ws->dof[k] = esmac_dof_fermions_alloc(par->n_sites, par->n_fermions);
      ws->rep[k] = esmac_rep_fermions_alloc(ws->dof[k]);
      ws->repinit[k] = esmac_rep_fermions_alloc(ws->dof[k]);
      esmac_multidx_set(ws->midx, k, esmac_dof_fermions_ns(ws->dof[k]));
      ws->ns = ws->ns * esmac_dof_fermions_ns(ws->dof[k]);
    }
    ws->maxrowlength = 2 * par->n_species * par->n_sites;
    if (info) {
      *info = 0;
    }
    return ws;
  } else {
    if (info) {
    *info = ESMAC_ERANGE;
    }
    return NULL;
  }
}

void esmac_matrix_FreeFermionChain_free(esmac_matrix_FreeFermionChain_work_t *ws) {
  if (ws) {
    if (ws->dof) {
      int k;
      for (k=0;k<ws->ndof;k++) {
        if (ws->dof[k]) {
          esmac_dof_fermions_free(ws->dof[k]);
        }
      }
      free(ws->dof);
    }
    if (ws->rep) {
      int k;
      for (k=0;k<ws->ndof;k++) {
        if (ws->rep[k]) {
          esmac_rep_fermions_free(ws->rep[k]);
        }
      }
      free(ws->rep);
    }
    if (ws->repinit) {
      int k;
      for (k=0;k<ws->ndof;k++) {
        if (ws->rep[k]) {
          esmac_rep_fermions_free(ws->repinit[k]);
        }
      }
      free(ws->repinit);
    }
    if (ws->midx) {esmac_multidx_free(ws->midx);}
    if (ws->row) {esmac_sparserow_real_free(ws->row);}
    free(ws);
  }
}

int esmac_matrix_FreeFermionChain_set_info(const esmac_matrix_FreeFermionChain_params_t * par, esmac_matrix_FreeFermionChain_work_t *ws, esmac_matrix_info_t *info) {
  info->nrow = ws->ns;
  info->ncol = ws->ns;
  info->maxrowlength = ws->maxrowlength;
  info->valtype = ESMAC_VAL_REAL;
  info->symmetry = ESMAC_SYMMETRIC;
  return 0;
}

int esmac_matrix_FreeFermionChain_row(const esmac_matrix_FreeFermionChain_params_t * par, esmac_matrix_FreeFermionChain_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val) {
  if (0 <= irow && irow < ws->ns) {
    if (par->t != 0.0) {
    
      int k;
      for (k=0;k<ws->ndof;k++) {
        esmac_fermions_decode(ws->dof[k], esmac_multidx_decode(ws->midx, k, irow), ws->repinit[k]);
      }
 
      esmac_sparserow_real_zero(ws->row);
       
      for (k=0;k<ws->ndof;k++) {
        int i;
        // hop
        for (i=0;i<par->n_sites-1;i++) {
          esmac_rep_fermions_copy(ws->dof[k], ws->repinit[k], ws->rep[k]);
          double val = esmac_op_fermions_hop(ws->dof[k], ws->rep[k], i,i+1);
          if (val != 0.0) {
            esmac_sparserow_real_add(ws->row, par->t * val, 
                                     esmac_multidx_upd(ws->midx, k, esmac_fermions_encode(ws->dof[k], ws->rep[k]),
                                                       irow));
          }
        }
        if (par->boundary_conditions == ESMAC_PBC) {
          esmac_rep_fermions_copy(ws->dof[k], ws->repinit[k], ws->rep[k]);
          double val = esmac_op_fermions_hop(ws->dof[k], ws->rep[k], 0,par->n_sites-1);
          if (val != 0.0) {
            esmac_sparserow_real_add(ws->row, par->t * val, 
                                     esmac_multidx_upd(ws->midx, k, esmac_fermions_encode(ws->dof[k], ws->rep[k]),
                                                       irow));
          }
        }
      }
     
 
      int nr = esmac_sparserow_real_to_idxval(ws->row, ws->maxrowlength, 0, cind, val);      
    
      return nr;
    } else {
      return 0;
    }
  } else {
    return ESMAC_ERANGE;
  }   
           
}
