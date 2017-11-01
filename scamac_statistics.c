#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>

#include "scamac_include.h"
#include "scamac_internal.h"
#include "scamac_aux.h"
#include "scamac_string.h"
#include "scamac_statistics.h"

ScamacErrorCode scamac_collect_statistics_and_pattern(const ScamacGenerator * gen, ScamacFlag flag, scamac_matrix_statistics_st * st, scamac_matrix_pattern_st ** pt) {
  ScamacErrorCode err;
  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EINVALID;
  }
  bool fl_do_stat = st != NULL;
  bool fl_do_pat  = pt != NULL;
  if (!(fl_do_stat || fl_do_pat)) {
    return SCAMAC_EOK;  // nothing to do
  }

  if (flag != SCAMAC_DEFAULT) {
    return SCAMAC_EFAIL;  // not yet implemented
  }

  ScamacIdx nrow = scamac_generator_query_nrow   (gen);
  ScamacIdx ncol = scamac_generator_query_ncol   (gen);
  int valtype    = scamac_generator_query_valtype(gen);

  if (fl_do_stat) {
    scamac_statistics_empty(st, nrow,ncol,valtype);
  }
  if (fl_do_pat) {
    err = scamac_pattern_alloc(1024,1024,pt);
    if (err) {
      return err|SCAMAC_EINTERNAL;
    }
    err = scamac_pattern_empty(*pt, nrow,ncol,valtype);
    if (err) {
      return err|SCAMAC_EINTERNAL;
    }
  }

  ScamacWorkspace * ws;

  err = scamac_workspace_alloc(gen, &ws);
  if (err) {
    return err|SCAMAC_EINTERNAL;
  }


  ScamacIdx *cind;
  double *val;
  err = scamac_alloc_cind_val(gen, flag, &cind, &val);
  if (err) {
    return err|SCAMAC_EINTERNAL;
  }

  ScamacIdx irow, nzr;
  for (irow=0; irow<nrow; irow++) {
    err = scamac_generate_row(gen, ws, irow, flag, &nzr, cind, val);
    if (err) {
      return err|SCAMAC_EINTERNAL;
    }
    if (fl_do_stat) {
      err = scamac_statistics_update(st, irow, nzr, cind, val);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
    }
    if (fl_do_pat) {
      err = scamac_pattern_update(*pt, irow, nzr, cind);
      if (err) {
        return err|SCAMAC_EINTERNAL;
      }
    }
  }

  free(cind);
  free(val);
  err = scamac_workspace_free(ws);
  if (err) {
    return err|SCAMAC_EINTERNAL;
  }

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_statistics_empty  (scamac_matrix_statistics_st * st, ScamacIdx nrow, ScamacIdx ncol, int valtype) {
  if (!st) {
    return SCAMAC_ENULL;
  }
  if ((nrow<1) || (ncol<1)) {
    return SCAMAC_ERANGE;
  }
  if ((valtype != SCAMAC_VAL_REAL) && (valtype != SCAMAC_VAL_COMPLEX)) {
    return SCAMAC_EUNKNOWN;
  }
  st->nrow=nrow;
  st->ncol=ncol;
  st->valtype=valtype;
  st->ncontributed=0;

  st->n_nz = 0;
  st->n_nz_left=0;
  st->n_nz_right=0;
  st->n_nz_row_min=0;
  st->n_nz_row_max=0;
  st->n_nz_row_max_left=0;
  st->n_nz_row_max_right=0;
  st->n_zero_row=0;
  st->n_zero_diag=0;
  st->bw_left=0;
  st->bw_right=0;
  st->bw=0;

  st->v_min_re=0.0;
  st->v_max_re=0.0;
  st->v_min_im=0.0;
  st->v_max_im=0.0;
  st->v_abs=0.0;
  st->v_min_re_diag=0.0;
  st->v_max_re_diag=0.0;
  st->v_min_im_diag=0.0;
  st->v_max_im_diag=0.0;
  st->v_abs_diag=0.0;
  st->n_diag_dominant=0;
  st->diag_minus_offdiag=0.0;

  // spectral statistics
  st->gershgorin_min_re=0.0;
  st->gershgorin_max_re=0.0;
  st->gershgorin_min_im=0.0;
  st->gershgorin_max_im=0.0;

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_statistics_update (scamac_matrix_statistics_st * st, ScamacIdx irow, ScamacIdx nzr, const ScamacIdx * cind, const double * val) {
  if (!st || !cind) {
    return SCAMAC_ENULL;
  }
  if (val) {
    if ((st->valtype != SCAMAC_VAL_REAL) && (st->valtype != SCAMAC_VAL_COMPLEX)) {
      return SCAMAC_EFAIL|SCAMAC_EINTERNAL;
    }
  }
  if ((irow<0) || (irow>=st->nrow)) {
    return SCAMAC_EINVALID;
  }

  ScamacIdx zero_diag=1;
  ScamacIdx row_left=0, row_right=0;
  double vd_re=0.0, vd_im=0.0, vsum=0.0;  // note: the Gershgorin estimate is not correct in case of doublets!
  ScamacIdx i;
  for (i=0; i<nzr; i++) {
    ScamacIdx c = cind[i];
    if (c<irow) {
      (st->n_nz_left)++;
      st->bw_left = MAX(st->bw_left,irow-c);
      row_left++;
    } else if (c>irow) {
      (st->n_nz_right)++;
      st->bw_right = MAX(st->bw_right,c-irow);
      row_right++;
    } else {
      zero_diag=0;
    }
    double v_re,v_im;
    if (st->valtype == SCAMAC_VAL_REAL) {
      v_re = val[i];
      v_im = 0.0;
    } else {// SCAMAC_VAL_COMPLEX
      v_re = val[2*i  ];
      v_im = val[2*i+1];
    }

    if (st->ncontributed==0) {
      st->v_min_re=v_re;
      st->v_max_re=v_re;
      st->v_min_im=v_im;
      st->v_max_im=v_im;
    } else {
      st->v_min_re=MIN(st->v_min_re,v_re);
      st->v_max_re=MAX(st->v_max_re,v_re);
      st->v_min_im=MIN(st->v_min_im,v_im);
      st->v_max_im=MAX(st->v_max_im,v_im);
    }
    if ((v_re != 0.0)||(v_im != 0.0)) {
      if (st->ncontributed==0) {
        st->v_abs=sqrt(v_re*v_re+v_im*v_im);
      } else {
        st->v_abs=MIN(st->v_abs,sqrt(v_re*v_re+v_im*v_im));
      }
    }
    if (c == irow) {
      if (st->ncontributed==0) {
        st->v_min_re_diag=v_re;
        st->v_max_re_diag=v_re;
        st->v_min_im_diag=v_im;
        st->v_max_im_diag=v_im;
      } else {
        st->v_min_re_diag=MIN(st->v_min_re_diag,v_re);
        st->v_max_re_diag=MAX(st->v_max_re_diag,v_re);
        st->v_min_im_diag=MIN(st->v_min_im_diag,v_im);
        st->v_max_im_diag=MAX(st->v_max_im_diag,v_im);
      }
      if ((v_re != 0.0) ||(v_im != 0.0)) {
        if (st->ncontributed==0) {
          st->v_abs_diag=sqrt(v_re*v_re+v_im*v_im);
        } else {
          st->v_abs_diag=MIN(st->v_abs,sqrt(v_re*v_re+v_im*v_im));
        }
      }
      vd_re = v_re;
      vd_im = v_im;
    } else {
      vsum = vsum + sqrt(v_re*v_re+v_im*v_im);
    }
  }
  st->n_nz = st->n_nz+nzr;
  if (st->ncontributed==0) {
    st->n_nz_row_min=nzr;
  } else {
    st->n_nz_row_min=MIN(st->n_nz_row_min,nzr);
  }
  st->n_nz_row_max=MAX(st->n_nz_row_max,nzr);
  st->n_nz_row_max_left =MAX(st->n_nz_row_max_left,row_left );
  st->n_nz_row_max_right=MAX(st->n_nz_row_max_right,row_right);
  if (nzr==0) {
    (st->n_zero_row)++ ;
  }
  if (zero_diag) {
    (st->n_zero_diag)++;
    /*
    st->vmin_diag=MIN(st->vmin_diag,0.0);
    st->vmax_diag=MAX(st->vmax_diag,0.0);
    st->veps_diag=0.0;
    */
  }
  if (st->ncontributed==0) {
    st->gershgorin_min_re=vd_re-vsum;
    st->gershgorin_max_re=vd_re+vsum;
    st->gershgorin_min_im=vd_im-vsum;
    st->gershgorin_max_im=vd_im+vsum;
  } else {
    st->gershgorin_min_re=MIN(st->gershgorin_min_re,vd_re-vsum);
    st->gershgorin_max_re=MAX(st->gershgorin_max_re,vd_re+vsum);
    st->gershgorin_min_im=MIN(st->gershgorin_min_im,vd_im-vsum);
    st->gershgorin_max_im=MAX(st->gershgorin_max_im,vd_im+vsum);
  }
  if (sqrt(vd_re*vd_re+vd_im*vd_im) > vsum) {
    st->n_diag_dominant++;
  }
  if (st->ncontributed==0) {
    st->diag_minus_offdiag = sqrt(vd_re*vd_re+vd_im*vd_im) - vsum;
  } else {
    st->diag_minus_offdiag = MIN(st->diag_minus_offdiag, sqrt(vd_re*vd_re+vd_im*vd_im) - vsum);
  }

  (st->ncontributed)++;

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_statistics_combine(scamac_matrix_statistics_st * stcomb, const scamac_matrix_statistics_st * st) {
  if (!st || !stcomb) {
    return SCAMAC_ENULL;
  }
  if ((st->nrow != stcomb->nrow) || (st->ncol != stcomb->ncol) || (st->valtype != stcomb->valtype)) {
    return SCAMAC_EINVALID;
  }

  if (st->ncontributed==0) {
    return SCAMAC_EOK;  // nothing to do
  }
  if (stcomb->ncontributed==0) {
    // simple copy will do
    memcpy(stcomb, st, sizeof *st);
    return SCAMAC_EOK;
  }

  stcomb->ncontributed += st->ncontributed;
  // don't check this .. it might be wrong for transposed of non-square matrix ... ???
  // if (stcomb->ncontributed > stcomb->nrow) { return SCAMAC_EFAIL; }

  stcomb->n_nz += st->n_nz;
  stcomb->n_nz_left += st->n_nz_left;
  stcomb->n_nz_right += st->n_nz_right;
  stcomb->n_nz_row_min=MIN(stcomb->n_nz_row_min,st->n_nz_row_min);
  stcomb->n_nz_row_max=MAX(stcomb->n_nz_row_max,st->n_nz_row_max);
  stcomb->n_nz_row_max_left=MAX(stcomb->n_nz_row_max_left,st->n_nz_row_max_left);
  stcomb->n_nz_row_max_right=MAX(stcomb->n_nz_row_max_right,st->n_nz_row_max_right);
  stcomb->n_zero_row += st->n_zero_row;
  stcomb->n_zero_diag += st->n_zero_diag;
  stcomb->bw_left=MAX(stcomb->bw_left,st->bw_left);
  stcomb->bw_right=MAX(stcomb->bw_right,st->bw_right);
  stcomb->bw=MAX(stcomb->bw_left,stcomb->bw_right);

  stcomb->v_min_re=MIN(stcomb->v_min_re,st->v_min_re);
  stcomb->v_max_re=MAX(stcomb->v_max_re,st->v_max_re);
  stcomb->v_min_im=MIN(stcomb->v_min_im,st->v_min_im);
  stcomb->v_max_im=MAX(stcomb->v_max_im,st->v_max_im);
  stcomb->v_abs=MAX(stcomb->v_abs,st->v_abs);
  stcomb->v_min_re_diag=MIN(stcomb->v_min_re_diag,st->v_min_re_diag);
  stcomb->v_max_re_diag=MAX(stcomb->v_max_re_diag,st->v_max_re_diag);
  stcomb->v_min_im_diag=MIN(stcomb->v_min_im_diag,st->v_min_im_diag);
  stcomb->v_max_im_diag=MAX(stcomb->v_max_im_diag,st->v_max_im_diag);
  stcomb->v_abs_diag=MAX(stcomb->v_abs_diag,st->v_abs_diag);
  stcomb->diag_minus_offdiag=MIN(stcomb->diag_minus_offdiag,st->diag_minus_offdiag);

  // spectral statistics
  stcomb->gershgorin_min_re=MIN(stcomb->gershgorin_min_re,st->gershgorin_min_re);
  stcomb->gershgorin_max_re=MAX(stcomb->gershgorin_max_re,st->gershgorin_max_re);
  stcomb->gershgorin_min_im=MIN(stcomb->gershgorin_min_im,st->gershgorin_min_im);
  stcomb->gershgorin_max_im=MAX(stcomb->gershgorin_max_im,st->gershgorin_max_im);

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_statistics_print  (const scamac_matrix_statistics_st * st, char ** desc) {
  if (!st || !desc) {
    return SCAMAC_ENULL;
  }
  const int nbuf = 256;
  char *buf = malloc( nbuf * sizeof *buf);

  scamac_string_st mystr;
  scamac_string_empty(&mystr);

  snprintf(buf,nbuf," === matrix ===\n");
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"rows:    %"SCAMACPRIDX"\n",st->nrow);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"columns: %"SCAMACPRIDX"\n",st->ncol);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"\n === pattern ===\n");
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements: %"SCAMACPRIDX"\n",st->n_nz);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements (left): %"SCAMACPRIDX"\n",st->n_nz_left);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements (right): %"SCAMACPRIDX"\n",st->n_nz_right);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements per row (min): %"SCAMACPRIDX"\n",st->n_nz_row_min);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements per row (max): %"SCAMACPRIDX"\n",st->n_nz_row_max);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements per row (max, left): %"SCAMACPRIDX"\n",st->n_nz_row_max_left);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"non-zero elements per row (max, right): %"SCAMACPRIDX"\n",st->n_nz_row_max_right);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"zero rows: %"SCAMACPRIDX"\n",st->n_zero_row);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"zero diagonal elements: %"SCAMACPRIDX"\n",st->n_zero_diag);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"bandwidth: %"SCAMACPRIDX"\n",st->bw);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"bandwidth (left): %"SCAMACPRIDX"\n",st->bw_left);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"bandwidth (right): %"SCAMACPRIDX"\n",st->bw_right);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"\n === values ===\n");
  scamac_string_append(&mystr, buf);
  if (st->valtype==SCAMAC_VAL_REAL) {
    snprintf(buf,nbuf,"minimal: %f\n",st->v_min_re);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"maximal: %f\n",st->v_max_re);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal |non-zero|: %f\n",st->v_abs);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal diagonal entry: %f\n",st->v_min_re_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"maximal diagonal entry: %f\n",st->v_max_re_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal |non-zero diagonal|: %f\n",st->v_abs_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"number of diagonally dominant rows: %"SCAMACPRIDX"\n",st->n_diag_dominant);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal |diag|-sum |offdiag|: %f\n",st->diag_minus_offdiag);
    scamac_string_append(&mystr, buf);
  } else if (st->valtype==SCAMAC_VAL_COMPLEX) {
    snprintf(buf,nbuf,"minimal real: %f\n",st->v_min_re);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"maximal real: %f\n",st->v_max_re);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal imag: %f\n",st->v_min_im);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"maximal imag: %f\n",st->v_max_im);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal |non-zero|: %f\n",st->v_abs);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal real diagonal entry: %f\n",st->v_min_re_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"maximal real diagonal entry: %f\n",st->v_max_re_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal imag diagonal entry: %f\n",st->v_min_im_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"maximal imag diagonal entry: %f\n",st->v_max_im_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal |non-zero diagonal|: %f\n",st->v_abs_diag);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"number of diagonally dominant rows: %"SCAMACPRIDX"\n",st->n_diag_dominant);
    scamac_string_append(&mystr, buf);
    snprintf(buf,nbuf,"minimal |diag|-sum |offdiag|: %f\n",st->diag_minus_offdiag);
    scamac_string_append(&mystr, buf);
  } else {
    return SCAMAC_EFAIL;
  }
  snprintf(buf,nbuf,"\n === spectrum ===\n");
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"Gershgorin estimate (minimal real): %f\n",st->gershgorin_min_re);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"Gershgorin estimate (maximal real): %f\n",st->gershgorin_max_re);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"Gershgorin estimate (minimal imag): %f\n",st->gershgorin_min_im);
  scamac_string_append(&mystr, buf);
  snprintf(buf,nbuf,"Gershgorin estimate (maximal imag): %f\n",st->gershgorin_max_im);
  scamac_string_append(&mystr, buf);

  *desc = scamac_string_get(&mystr);

  return SCAMAC_EOK;
}


ScamacErrorCode scamac_pattern_alloc  (int px, int py, scamac_matrix_pattern_st ** pt) {
  if ((px<1) || (py<1)) {
    return SCAMAC_ERANGE;
  }
  if ((px>1024) || (py>1024)) {
    return SCAMAC_EWARNING;  // rather large, isn't it?
  }
  if (!pt) {
    return SCAMAC_ENULL;
  }
  scamac_matrix_pattern_st *my_pt = malloc(sizeof *my_pt);
  if (!my_pt) {
    return SCAMAC_EMALLOCFAIL;
  }
  my_pt->px = px;
  my_pt->py = py;
  my_pt->pat = malloc(px*py * sizeof *(my_pt->pat));
  if (!my_pt->pat) {
    return SCAMAC_EMALLOCFAIL;
  }
  my_pt->nrow=0;
  my_pt->ncol=0;
  *pt = my_pt;
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_pattern_empty  (scamac_matrix_pattern_st * pt, ScamacIdx nrow, ScamacIdx ncol, int valtype) {
  if (!pt) {
    return SCAMAC_ENULL;
  }
  if (!pt->pat) {
    return SCAMAC_ENULL;
  }
  if ((nrow<1) || (ncol<1)) {
    return SCAMAC_ERANGE;
  }
  pt->nrow=nrow;
  pt->ncol=ncol;
  pt->ncontributed=0;
  // memset
  memset(pt->pat,0,pt->px * pt->py * sizeof *(pt->pat));
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_pattern_update (scamac_matrix_pattern_st * pt, ScamacIdx irow, ScamacIdx nzr, const ScamacIdx * cind) {
  if (!pt || !cind) {
    return SCAMAC_ENULL;
  }
  /*
  if (pt->nrow == 0) {
    pt->nrow = nrow;
  } else {
    if (pt->nrow != nrow) { return SCAMAC_EINVALID; }
  }
  if (pt->ncol == 0) {
    pt->ncol = ncol;
  } else {
    if (pt->ncol != ncol) { return SCAMAC_EINVALID; }
  }
  */
  if ((irow<0) || (irow>=pt->nrow)) {
    return SCAMAC_EINVALID;
  }

  int y = floor (  (double) irow * (double) pt->py / (double) (pt->nrow) );
  ScamacIdx i;
  for (i=0; i<nzr; i++) {
    if ((cind[i]<0) || (cind[i]>=pt->ncol)) {
      return SCAMAC_EFAIL;
    }
    int x = floor (  (double) cind[i] * (double) pt->px / (double) (pt->ncol) );
    pt->pat[y*pt->px + x]=1;
    //    pt->pat[y*pt->px + x]++;
  }

  (pt->ncontributed)++;

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_pattern_combine(scamac_matrix_pattern_st * ptcomb, const scamac_matrix_pattern_st * pt) {
  if (!ptcomb || !pt) {
    return SCAMAC_ENULL;
  }
  if ((pt->px != ptcomb->px) || (pt->py != ptcomb->py)) {
    return SCAMAC_EINVALID;
  }
  if ((pt->nrow != ptcomb->nrow) || (pt->ncol != ptcomb->ncol)) {
    return SCAMAC_EINVALID;
  }

  ptcomb->ncontributed += pt->ncontributed;
  ScamacIdx i;
  for (i=0; i<pt->px*pt->py; i++) {
    ptcomb->pat[i] += pt->pat[i];
  }

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_pattern_print  (const scamac_matrix_pattern_st * pt, char ** desc) {
  if (!pt || !desc) {
    return SCAMAC_ENULL;
  }
  // 64*64 characters for pattern (with border: 66*66 for the "picture", to be precise, with newline: 67*66)
  char * pic = malloc(67*67* sizeof *pic);
  if (!pic) {
    return SCAMAC_EMALLOCFAIL;
  }
  int ix,iy;

  // border
  for (ix=1; ix<65; ix++) {
    pic[67*0 + ix]='-';
    pic[67*65 + ix]='-';
  }
  for (iy=1; iy<65; iy++) {
    pic[67*iy + 0]='|';
    pic[67*iy + 65]='|';
  }
  pic[67*0 + 0]='+';
  pic[67*65 + 0]='+';
  pic[67*0 + 65]='+';
  pic[67*65 + 65]='+';
  // newlinw
  for (iy=0; iy<66; iy++) {
    pic[67*iy + 66]='\n';
  }

  // pattern: initialize as blank
  for (iy=0; iy<64; iy++) {
    for (ix=0; ix<64; ix++) {
      pic[67*(iy+1)+(ix+1)]=' ';
    }
  }
  // set "dots"
  for (iy=0; iy<pt->py; iy++) {
    int y = floor( (double) iy * ( (double) 64 / (double) (pt->py) ) );
    for (ix=0; ix<pt->px; ix++) {
      int x = floor( (double) ix * ( (double) 64 / (double) (pt->px) ) );
      if (pt->pat[iy*pt->py+ix]) {
        // * for non-zero
        // o for non-zero on main diagonal (note: due reduced resolution, the diagonal "o" entries represents *a block of entries* of the sparsity pattern)
        // if (x==y) {
        if (ix==iy) {
          pic[67*(y+1)+(x+1)] = 'o';
        } else {
          pic[67*(y+1)+(x+1)] = '*';
        }
      }
    }
  }

  *desc = pic;
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_pattern_free   (scamac_matrix_pattern_st * pt) {
  if (pt) {
    free(pt->pat);
    free(pt);
  }
  return SCAMAC_EOK;
}
