#include <stdlib.h>

#include "esmac_constants.h"
#include "esmac_matrix_SpinChainXXZ.h"


esmac_matrix_SpinChainXXZ_work_t * esmac_matrix_SpinChainXXZ_alloc(const esmac_matrix_SpinChainXXZ_params_t * par, int * info) {
  esmac_matrix_SpinChainXXZ_work_t * ws = malloc(sizeof *ws);
  ws->row = esmac_sparserow_real_alloc();
  ws->dof = esmac_dof_spinssz_alloc(par->n_sites, par->n_up);
  ws->rep = esmac_rep_spins_alloc(ws->dof);
  ws->repinit = esmac_rep_spins_alloc(ws->dof);
  ws->ns = esmac_dof_spins_ns(ws->dof);
  ws->maxrowlength = par->n_sites + 1;
  if (info) {
    *info = 0;
  }
  return ws;
}

void esmac_matrix_SpinChainXXZ_free(esmac_matrix_SpinChainXXZ_work_t *ws) {
  if (ws) {
    if (ws->dof) {
      esmac_dof_spins_free(ws->dof);
    }
    if (ws->rep) {
      esmac_rep_spins_free(ws->rep);
    }
    if (ws->repinit) {
      esmac_rep_spins_free(ws->repinit);
    }
    if (ws->row) {esmac_sparserow_real_free(ws->row);}
    free(ws);
  }
}

int esmac_matrix_SpinChainXXZ_set_info(const esmac_matrix_SpinChainXXZ_params_t * par, esmac_matrix_SpinChainXXZ_work_t *ws, esmac_matrix_info_t *info) {
  info->nrow = ws->ns;
  info->ncol = ws->ns;
  info->maxrowlength = ws->maxrowlength;
  info->valtype = ESMAC_VAL_REAL;
  info->symmetry = ESMAC_SYMMETRIC;
  return 0;
}

int esmac_matrix_SpinChainXXZ_row(const esmac_matrix_SpinChainXXZ_params_t * par, esmac_matrix_SpinChainXXZ_work_t *ws,
				 esmac_idx_t irow, esmac_idx_t *cind, double *val) {
  if (0 <= irow && irow < ws->ns) {
       
    esmac_spins_decode(ws->dof, irow, ws->repinit);
     
    esmac_sparserow_real_zero(ws->row);
      
    if (par->Jz != 0.0) {
      int i;
      double val=0.0;
      for (i=0;i<par->n_sites-1;i++) {
        val += esmac_op_spins_SzSz(ws->dof, ws->repinit, i,i+1);
      }
      if (par->boundary_conditions == ESMAC_PBC) {
        val += esmac_op_spins_SzSz(ws->dof, ws->repinit, 0,par->n_sites-1);
      }
      esmac_sparserow_real_add(ws->row, par->Jz * val, irow);
    }
    
    if (par->Bz != 0.0) {
      int i;
      double val=0.0;
      for (i=0;i<par->n_sites;i++) {
        val += esmac_op_spins_Sz(ws->dof, ws->repinit, i);
      }
      esmac_sparserow_real_add(ws->row, par->Bz * val, irow);
    }
    
    if (par->Jxy != 0.0) {
      int i;
      for (i=0;i<par->n_sites-1;i++) {
        esmac_rep_spins_copy(ws->dof, ws->repinit, ws->rep);
        double val = esmac_op_spins_Sflip(ws->dof, ws->rep, i,i+1);
        if (val != 0.0) {
          esmac_sparserow_real_add(ws->row, par->Jxy * val, esmac_spins_encode(ws->dof, ws->rep));
        }
      }
      if (par->boundary_conditions == ESMAC_PBC) {
        esmac_rep_spins_copy(ws->dof, ws->repinit, ws->rep);
        double val = esmac_op_spins_Sflip(ws->dof, ws->rep, 0,par->n_sites-1);
        if (val != 0.0) {
          esmac_sparserow_real_add(ws->row, par->Jxy * val, esmac_spins_encode(ws->dof, ws->rep));
        }
      }
    }
    
    int nr = esmac_sparserow_real_to_idxval(ws->row, ws->maxrowlength, 0, cind, val);      
    
    return nr;
  
  } else {
    return ESMAC_ERANGE;
  }   
           
}
