#include <stdlib.h>

#include "esmac_constants.h"
#include "esmac_dof_fermions.h"
#include "esmac_mstate.h"
#include "esmac_matrix_Hubbard.h"


esmac_matrix_Hubbard_work_t * esmac_matrix_Hubbard_alloc(const esmac_matrix_Hubbard_params_t * par, int * info) {
  esmac_matrix_Hubbard_work_t * ws = malloc(sizeof *ws);
  ws->midx = esmac_multidx_alloc(2);
  ws->row = esmac_sparserow_real_alloc();
  ws->dof = malloc(2 * sizeof *(ws->dof));
  ws->rep = malloc(2 * sizeof *(ws->rep));
  ws->repinit = malloc(2 * sizeof *(ws->repinit));
  ws->ns = 1;
  int k;
  for (k=0;k<2;k++) {
    ws->dof[k] = esmac_dof_fermions_alloc(par->n_sites, par->n_fermions);
    ws->rep[k] = esmac_rep_fermions_alloc(ws->dof[k]);
    ws->repinit[k] = esmac_rep_fermions_alloc(ws->dof[k]);
    esmac_multidx_set(ws->midx, k, esmac_dof_fermions_ns(ws->dof[k]));
    ws->ns = ws->ns * esmac_dof_fermions_ns(ws->dof[k]);
  }
  ws->maxrowlength = 2 * par->n_sites + 1;
  ws->rg=NULL;
  ws->onsite=NULL;
  if (par->ranpot > 0.0) {
    ws->rg = esmac_rng_alloc(par->seed);
    ws->onsite = malloc(par->n_sites * sizeof *(ws->onsite));
    int i;
    for (i=0;i<par->n_sites;i++) {
      ws->onsite[i] = esmac_rng_get_double(ws->rg, -par->ranpot, par->ranpot, i);
    }
  }
  if (info) {
    *info = 0;
  }
  return ws;
}

void esmac_matrix_Hubbard_free(esmac_matrix_Hubbard_work_t *ws) {
  if (ws) {
    if (ws->dof) {
      int k;
      for (k=0;k<2;k++) {
        if (ws->dof[k]) {
          esmac_dof_fermions_free(ws->dof[k]);
        }
      }
      free(ws->dof);
    }
    if (ws->rep) {
      int k;
      for (k=0;k<2;k++) {
        if (ws->rep[k]) {
          esmac_rep_fermions_free(ws->rep[k]);
        }
      }
      free(ws->rep);
    }
    if (ws->repinit) {
      int k;
      for (k=0;k<2;k++) {
        if (ws->rep[k]) {
          esmac_rep_fermions_free(ws->repinit[k]);
        }
      }
      free(ws->repinit);
    }
    if (ws->midx) {esmac_multidx_free(ws->midx);}
    if (ws->row) {esmac_sparserow_real_free(ws->row);}
    if (ws->onsite) {free(ws->onsite);}
    if (ws->rg) {esmac_rng_free(ws->rg);}
    free(ws);
  }
}

int esmac_matrix_Hubbard_set_info(const esmac_matrix_Hubbard_params_t * par, esmac_matrix_Hubbard_work_t *ws, esmac_matrix_info_t *info) {
  info->nrow = ws->ns;
  info->ncol = ws->ns;
  info->maxrowlength = ws->maxrowlength;
  info->valtype = ESMAC_VAL_REAL;
  info->symmetry = ESMAC_SYMMETRIC;
  return 0;
}

int esmac_matrix_Hubbard_row(const esmac_matrix_Hubbard_params_t * par, esmac_matrix_Hubbard_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val) {
  if (0 <= irow && irow < ws->ns) {
       
    int k;
    for (k=0;k<2;k++) {
      esmac_fermions_decode(ws->dof[k], esmac_multidx_decode(ws->midx, k, irow), ws->repinit[k]);
    }
 
    esmac_sparserow_real_zero(ws->row);
      
    if (par->t != 0.0) {
       
      for (k=0;k<2;k++) {
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
      
    }
     
     // diagonal Hubbard term
    if (par->U != 0.0) {
      double val = 0.0;
      int i;
      for (i=0;i<par->n_sites;i++) {
        val = val + esmac_op_fermions_cdc(ws->dof[0], ws->repinit[0], i) * esmac_op_fermions_cdc(ws->dof[1], ws->repinit[1], i);
      }
      esmac_sparserow_real_add(ws->row, par->U * val, irow);
    }
    
    // diagonal on-site potential
    if (ws->onsite) {
      double val = 0.0;
      int i;
      for (i=0;i<par->n_sites;i++) {
        val = val + ws->onsite[i]*(esmac_op_fermions_cdc(ws->dof[0], ws->repinit[0], i) 
                                   + esmac_op_fermions_cdc(ws->dof[1], ws->repinit[1], i));
      }
      esmac_sparserow_real_add(ws->row, val, irow);
    } 
    
    int nr = esmac_sparserow_real_to_idxval(ws->row, ws->maxrowlength, 0, cind, val);      
    
    return nr;
  
  } else {
    return ESMAC_ERANGE;
  }   
           
}
