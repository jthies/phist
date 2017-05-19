#include <stdlib.h>

#include "esmac_constants.h"
#include "esmac_matrix_FreeBosonChain.h"


esmac_matrix_FreeBosonChain_work_t * esmac_matrix_FreeBosonChain_alloc(const esmac_matrix_FreeBosonChain_params_t * par, int * info) {
  esmac_matrix_FreeBosonChain_work_t * ws = malloc(sizeof *ws);
  if (par->n_species != 0) {
    ws->ndof = par->n_species;
    ws->dof = malloc(par->n_species * sizeof *(ws->dof));
    ws->ms = esmac_mstate_alloc(par->n_species);
    ws->ns = 1;
    int k;
    for (k=0;k<par->n_species;k++) {
      ws->dof[k] = esmac_dof_bosons_alloc(par->n_sites, par->n_bosons, par->n_bosons);
      //esmac_mstate_set_qs(ws->ms, k,  esmac_dof_fermions_ns(ws->dof[k]), esmac_dof_fermions_qs(ws->dof[k]));
      esmac_mstate_register_dof_bosons(ws->ms, k,  ws->dof[k]);
      ws->ns = ws->ns * esmac_dof_bosons_ns(ws->dof[k]);
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

void esmac_matrix_FreeBosonChain_free(esmac_matrix_FreeBosonChain_work_t *ws) {
  if (ws) {
    if (ws->dof) {
      int k;
      for (k=0;k<ws->ndof;k++) {
        if (ws->dof[k]) {
          esmac_dof_bosons_free(ws->dof[k]);
        }
      }
      free(ws->dof);
    }
    if (ws->ms) {esmac_mstate_free(ws->ms);}
    free(ws);
  }
}

int esmac_matrix_FreeBosonChain_set_info(const esmac_matrix_FreeBosonChain_params_t * par, esmac_matrix_FreeBosonChain_work_t *ws, esmac_matrix_info_t *info) {
  info->nrow = ws->ns;
  info->ncol = ws->ns;
  info->maxrowlength = ws->maxrowlength;
  info->valtype = ESMAC_VAL_REAL;
  info->symmetry = ESMAC_SYMMETRIC;
  return 0;
}

int esmac_matrix_FreeBosonChain_row(const esmac_matrix_FreeBosonChain_params_t * par, esmac_matrix_FreeBosonChain_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val) {
  if (0 <= irow && irow < ws->ns) {
    esmac_mstate_decode(ws->ms, irow);
    if (par->t != 0.0) {
      int k;
      for (k=0;k<ws->ndof;k++) {
        int i;
        // hop
        for (i=0;i<par->n_sites-1;i++) {
          esmac_op_bosons_bdibj(ws->dof[k],  par->t, i,i+1);
          esmac_op_bosons_bdibj(ws->dof[k],  par->t, i+1,i);
        }
        if (par->boundary_conditions == ESMAC_PBC) {
          esmac_op_bosons_bdibj(ws->dof[k],  par->t, 0,par->n_sites-1);
          esmac_op_bosons_bdibj(ws->dof[k],  par->t, par->n_sites-1,0);
        }
        // esmac_mstate_row_add(ws->ms, 1.0);
        esmac_mstate_row_add_one_qs(ws->ms, k, 1.0);
      }
    }
 
    int nr = esmac_mstate_row_to_idxval(ws->ms, ws->maxrowlength, 0, cind, val);      
    
    return nr;

  } else {
    return ESMAC_ERANGE;
  }   
           
}
