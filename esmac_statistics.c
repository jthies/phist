#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "esmac_types.h"
#include "esmac_statistics.h"

void esmac_print_statistics(FILE *f, const esmac_matrix_statistics_t *st) {
  fprintf(f," === matrix ===\n");
  fprintf(f,"rows:    %"ESMACPRIDX"\n",st->nr);
  fprintf(f,"columns: %"ESMACPRIDX"\n",st->nc);
  fprintf(f,"\n === pattern ===\n");
  fprintf(f,"non-zero elements: %"ESMACPRIDX"\n",st->n_nz);
  fprintf(f,"non-zero elements (left): %"ESMACPRIDX"\n",st->n_nz_left);
  fprintf(f,"non-zero elements (right): %"ESMACPRIDX"\n",st->n_nz_right);
  fprintf(f,"non-zero elements per row (min): %"ESMACPRIDX"\n",st->n_nz_row_min);
  fprintf(f,"non-zero elements per row (max): %"ESMACPRIDX"\n",st->n_nz_row_max);
  fprintf(f,"non-zero elements per row (max, left): %"ESMACPRIDX"\n",st->n_nz_row_max_left);
  fprintf(f,"non-zero elements per row (max, right): %"ESMACPRIDX"\n",st->n_nz_row_max_right);
  fprintf(f,"zero rows: %"ESMACPRIDX"\n",st->n_zero_row);
  fprintf(f,"zero diagonal elements: %"ESMACPRIDX"\n",st->n_zero_diag);
  fprintf(f,"bandwidth: %"ESMACPRIDX"\n",st->bw);
  fprintf(f,"bandwidth (left): %"ESMACPRIDX"\n",st->bw_left);
  fprintf(f,"bandwidth (right): %"ESMACPRIDX"\n",st->bw_right);
  fprintf(f,"\n === values ===\n");
  fprintf(f,"minimal: %f\n",st->vmin);
  fprintf(f,"maximal: %f\n",st->vmax);
  fprintf(f,"minimal |non-zero|: %f\n",st->veps);
  fprintf(f,"minimal diagonal entry: %f\n",st->vmin_diag);
  fprintf(f,"maximal diagonal entry: %f\n",st->vmax_diag);
  fprintf(f,"minimal |non-zero diagonal|: %f\n",st->veps_diag);
  fprintf(f,"\n === spectrum ===\n");
  fprintf(f,"Gershgorin estimate (minimal): %f\n",st->gershgorin_min);
  fprintf(f,"Gershgorin estimate (maximal): %f\n",st->gershgorin_max);
}


void esmac_collect_matrix_statistics(const esmac_generator_t * gen, esmac_matrix_statistics_t *st, int pattern_px, int *pattern, bool show_progress) {
  int info;
  esmac_generator_work_t * my_ws = esmac_generator_alloc(gen, &info);
 
  esmac_idx_t nrow = esmac_generator_query(gen,my_ws,"nrow");
  esmac_idx_t ncol = esmac_generator_query(gen,my_ws,"ncol"); 
  esmac_idx_t mrow = esmac_generator_query(gen,my_ws,"maxnzrow");

  esmac_idx_t *cind = malloc(mrow * sizeof *cind);
  double *val = malloc(mrow * sizeof *val);
  
  esmac_idx_t idx;
  
  st->nr = nrow;
  st->nc = ncol;
  
  st->n_nz = 0;
  st->n_nz_left=0;
  st->n_nz_right=0;
  st->n_nz_row_min=ncol+1; // something ridiculously large
  st->n_nz_row_max=0;
  st->n_nz_row_max_left=0;
  st->n_nz_row_max_right=0;
  st->n_zero_row=0;
  st->n_zero_diag=0;
  st->bw_left=0;
  st->bw_right=0;
  
  st->vmin=DBL_MAX;
  st->vmax=-DBL_MAX;
  st->veps=DBL_MAX;
  st->vmin_diag=DBL_MAX;
  st->vmax_diag=-DBL_MAX;
  st->veps_diag=DBL_MAX;

  // spectral statistics
  st->gershgorin_min=DBL_MAX;
  st->gershgorin_max=-DBL_MAX;
    
  for (idx=0;idx<nrow;idx++) {
    if (show_progress) {
      if ((idx > 0) && !(idx % 10000)) {
        fprintf(stdout,"[%s: generate row %"ESMACPRIDX" (%.2d%%)]\r",__func__,idx,(int) (100.0*idx/nrow));
        fflush(stdout);
      }
    }
    int k = esmac_generator_row(gen, my_ws, idx, cind, val);
    //beware about int <-> esmac_idx_t
    //  memcpy(&(sm->cind[n]), cind, k * sizeof *cind);
    int i;
    int zero_diag=1;
    int row_left=0, row_right=0;
    double vd=0.0, vsum=0.0;
    int y = floor (  (double) idx * (double) pattern_px / (double) (nrow) );
    // note: the Gershgorin estimate is not correct in case of doublets!
    for (i=0;i<k;i++) {
      int c = cind[i];
      int x = floor (  (double) cind[i] * (double) pattern_px / (double) (ncol) );
      if (pattern) {pattern[y*pattern_px+x]=1;}
      if (c<idx) { (st->n_nz_left)++; st->bw_left = MAX(st->bw_left,idx-c); row_left++; }
      else if (c>idx) { (st->n_nz_right)++; st->bw_right = MAX(st->bw_right,c-idx); row_right++; }
      else {zero_diag=0;}
      double v = val[i];
      st->vmin=MIN(st->vmin,v);
      st->vmax=MAX(st->vmax,v);
      if (v != 0.0) {st->veps=MIN(st->veps,fabs(v));}
      if (c == idx) {
        st->vmin_diag=MIN(st->vmin_diag,v);
        st->vmax_diag=MAX(st->vmax_diag,v);
        if (v != 0.0) {st->veps_diag=MIN(st->veps_diag,fabs(v));}
        vd = v;
      } else {
        vsum = vsum + fabs(v);
      }
    }
    st->n_nz = st->n_nz+k;
    st->n_nz_row_min=MIN(st->n_nz_row_min,k);
    st->n_nz_row_max=MAX(st->n_nz_row_max,k);
    st->n_nz_row_max_left =MAX(st->n_nz_row_max_left ,row_left );
    st->n_nz_row_max_right=MAX(st->n_nz_row_max_right,row_right);
    if (k==0) { (st->n_zero_row)++ ; }
    if (zero_diag) {
      (st->n_zero_diag)++;
      st->vmin_diag=MIN(st->vmin_diag,0.0);
      st->vmax_diag=MAX(st->vmax_diag,0.0);
      st->veps_diag=0.0;
    }
    st->gershgorin_min=MIN(st->gershgorin_min,vd-vsum);
    st->gershgorin_max=MAX(st->gershgorin_max,vd+vsum);
  }
  if (show_progress) {printf("\n\n");}
 
  st->bw = MAX(st->bw_left, st->bw_right);
  
  free(cind);
  free(val);
  esmac_generator_free(my_ws);
    
}
