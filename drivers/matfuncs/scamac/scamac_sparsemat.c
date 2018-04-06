#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "scamac_internal.h"
#include "scamac_aux.h"
#include "scamac_safeint.h"
#include "scamac_sparsemat.h"

ScamacErrorCode scamac_sparsemat_alloc(ScamacIdx nr, ScamacIdx nc, ScamacIdx ne, int valtype, scamac_sparsemat_st ** sm) {
  if ((valtype != SCAMAC_VAL_REAL) && (valtype != SCAMAC_VAL_COMPLEX)) {
    return SCAMAC_EINVALID;
  }
  if (!sm) {
    return SCAMAC_ENULL;
  }
  if ((nr<1) || (nc<1) || (ne<1)) {
    return SCAMAC_ERANGE;
  }

  scamac_sparsemat_st * my_sm = malloc(sizeof * my_sm);
  if (!my_sm) {
    return SCAMAC_EMALLOCFAIL;
  }
  my_sm->nr = nr;
  my_sm->nc = nc;
  my_sm->ne=0;
  my_sm->nemax = ne;
  my_sm->rptr = malloc((nr+1) * sizeof *(my_sm->rptr));
  if (!my_sm->rptr) {
    return SCAMAC_EMALLOCFAIL;
  }
  my_sm->cind = malloc( ne    * sizeof *(my_sm->cind));
  if (!my_sm->cind) {
    return SCAMAC_EMALLOCFAIL;
  }
  if (valtype == SCAMAC_VAL_REAL) {
    my_sm->val  = malloc(   ne  * sizeof *(my_sm->val) );
  } else if (valtype == SCAMAC_VAL_COMPLEX) {
    my_sm->val  = malloc( 2*ne  * sizeof *(my_sm->val) );
  } else {
    my_sm->val = NULL;
  }
  if (!my_sm->val) {
    return SCAMAC_EMALLOCFAIL;
  }
  my_sm->valtype = valtype;

  *sm = my_sm;
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_sparsemat_free(scamac_sparsemat_st * sm) {
  if (sm) {
    free(sm->rptr);
    free(sm->cind);
    free(sm->val );
    free(sm);
  }
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_sparsemat_from_generator(const ScamacGenerator * gen, scamac_sparsemat_st ** sm) {
  if (!gen || !sm) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }

  ScamacErrorCode err;
  ScamacWorkspace * my_ws;
  err = scamac_workspace_alloc(gen, &my_ws);
  if (err) {
    return err|SCAMAC_EINTERNAL;
  }

  ScamacIdx nrow = scamac_generator_query_nrow      (gen);
  ScamacIdx ncol = scamac_generator_query_ncol      (gen);
  ScamacIdx mrow = scamac_generator_query_maxnzrow  (gen);
  ScamacIdx valtype = scamac_generator_query_valtype(gen);
  if ( (nrow<1) || (ncol<1) || (mrow<1) ) {
    return SCAMAC_ERANGE|SCAMAC_EINTERNAL;
  }
  if ((valtype != SCAMAC_VAL_REAL) && (valtype != SCAMAC_VAL_COMPLEX)) {
    return SCAMAC_EINVALID|SCAMAC_EINTERNAL;
  }

  ScamacIdx *cind = malloc(mrow * sizeof *cind);
  if (!cind) {
    return SCAMAC_EMALLOCFAIL;
  }
  double *val = NULL;
  if (valtype == SCAMAC_VAL_REAL) {
    val = malloc(  mrow * sizeof *val);
  } else if (valtype == SCAMAC_VAL_COMPLEX) {
    val = malloc(2*mrow * sizeof *val);
  }
  if (!val) {
    return SCAMAC_EMALLOCFAIL;
  }

  ScamacIdx idx;

  // count # non-zeros
  ScamacIdx ne = 0;
  for (idx=0; idx<nrow; idx++) {
    ScamacIdx k;
    err = scamac_generate_row(gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val);
    if (err) {
      return err|SCAMAC_EINTERNAL;
    }
    ne = scamac_safe_add(ne,k); // ne=ne+k;
    if (ne<0) {
      return SCAMAC_EOVERFLOW;
    }
  }

  scamac_sparsemat_st * my_sm;
  err = scamac_sparsemat_alloc(nrow,ncol,ne,valtype,&my_sm);
  if (err) {
    return err|SCAMAC_EINTERNAL;
  }

  // create matrix
  ScamacIdx n = 0;
  my_sm->rptr[0]=0;
  for (idx=0; idx<nrow; idx++) {
    ScamacIdx k;
    err = scamac_generate_row(gen, my_ws, idx, SCAMAC_DEFAULT, &k, cind, val);
    if (err) {
      return err|SCAMAC_EINTERNAL;
    }
    //beware about int <-> ScamacIdx
    //  memcpy(&(sm->cind[n]), cind, k * sizeof *cind);
    int i;
    for (i=0; i<k; i++) {
      my_sm->cind[i+n]=cind[i];
    }
    if (valtype == SCAMAC_VAL_REAL) {
      memcpy(&(my_sm->val [n]), val, k * sizeof *val);
    } else if (valtype == SCAMAC_VAL_COMPLEX) {
      memcpy(&(my_sm->val [2*n]), val, 2*k * sizeof *val);
    }
    n=n+k;
    my_sm->rptr[idx+1]=n;
  }
  my_sm->ne=ne;

  free(cind);
  free(val);
  err = scamac_workspace_free(my_ws);
  if (err) {
    return err|SCAMAC_EINTERNAL;
  }

  *sm = my_sm;

  return SCAMAC_EOK;

}


ScamacErrorCode scamac_sparsemat_mvm(const scamac_sparsemat_st * sm, const double * x, double * y, double alpha, double beta, double gamma) {
  if (!sm || !x || !y) {
    return SCAMAC_ENULL;
  }
  if (sm->valtype != SCAMAC_VAL_REAL) {
    return SCAMAC_EINVALID;
  }
  if (gamma != 0.0 && sm->nr != sm->nc) {
    return SCAMAC_EINVALID;  // gamma != 0 requires SQUARE matrix a.nr == a.nc\n!
  }


  ScamacIdx i,j;

  // check for early return. TODO: call BLAS ?
  if (alpha==0.0 && beta==0.0) {
    for (i=0; i<sm->nr; i++) {
      y[i]=gamma*x[i];
    }
    return SCAMAC_EOK;
  }
  if (alpha==0.0) {
    for (i=0; i<sm->nr; i++) {
      y[i]=beta*y[i]+gamma*x[i];
    }
    return SCAMAC_EOK;
  }

  for (i=0; i<sm->nr; i++) {
    //if beta=0 but y[i]=NaN an error might occur because 0*NaN /= 0
    //therefore: check
    if (beta==0.0) {
      y[i]=gamma*x[i];
    } else {
      y[i]=beta*y[i]+gamma*x[i];
    }
    for (j=sm->rptr[i]; j<sm->rptr[i+1]; j++) {
      y[i]=y[i]+alpha*sm->val[j]*x[sm->cind[j]];
    }
  }

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_sparsemat_mvm_cplx(const scamac_sparsemat_st *sm, const double complex * x, double complex * y,
    double complex alpha, double complex beta, double complex gamma) {
  if (!sm || !x || !y) {
    return SCAMAC_ENULL;
  }
  if (sm->valtype != SCAMAC_VAL_COMPLEX) {
    return SCAMAC_EINVALID;
  }
  if (gamma != 0.0 && sm->nr != sm->nc) {
    return SCAMAC_EINVALID;  // gamma != 0 requires SQUARE matrix a.nr == a.nc\n!
  }

  ScamacIdx i,j;

  // check for early return. TODO: call BLAS ?
  if (alpha==0.0 && beta==0.0) {
    for (i=0; i<sm->nr; i++) {
      y[i]=gamma*x[i];
    }
    return SCAMAC_EOK;
  }
  if (alpha==0.0) {
    for (i=0; i<sm->nr; i++) {
      y[i]=beta*y[i]+gamma*x[i];
    }
    return SCAMAC_EOK;
  }

  double complex *myval = (double complex *) sm->val;

  for (i=0; i<sm->nr; i++) {
    //if beta=0 but y[i]=NaN an error might occur because 0*NaN /= 0
    //therefore: check
    if (beta==0.0) {
      y[i]=gamma*x[i];
    } else {
      y[i]=beta*y[i]+gamma*x[i];
    }
    for (j=sm->rptr[i]; j<sm->rptr[i+1]; j++) {
      y[i]=y[i]+alpha*myval[j]*x[sm->cind[j]];
    }
  }

  return SCAMAC_EOK;
}

ScamacIdx scamac_sparsemat_maxrowlength(const scamac_sparsemat_st * sm) {
  if (!sm) {
    return 0;
  } else {
    ScamacIdx mrwl = 0;
    ScamacIdx i;
    for (i=0; i<sm->nr; i++) {
      mrwl = MAX(mrwl,sm->rptr[i+1]-sm->rptr[i]);
    }
    return mrwl;
  }
}


ScamacErrorCode scamac_sparsemat_check_symmetry(const scamac_sparsemat_st * sm, bool * symm_pattern, bool * symm_value) {
  if (!sm) {
    return SCAMAC_ENULL;
  }
  if (sm->nr != sm->nc) {
    if (symm_pattern) {
      *symm_pattern=false;
    }
    if (symm_value) {
      *symm_value=false;
    }
    return SCAMAC_EINVALID; // check requires SQUARE matrix
  }

  if (!symm_pattern && !symm_value) {
    return SCAMAC_EOK;  // nothing to do here
  }

  bool my_s_pat = true;
  bool my_s_val = true;


  // ScamacIdx *ivec = calloc(sm->nr, sizeof *ivec);
  ScamacIdx *ivec = malloc(sm->nr * sizeof *ivec);
  if (!ivec) {
    return SCAMAC_EMALLOCFAIL;
  }
  memcpy(ivec, sm->rptr, sm->nr * sizeof *ivec);

  ScamacIdx i,j;
  for (i=0; i<sm->nr; i++) {
    for (j=sm->rptr[i]; j<sm->rptr[i+1]; j++) {
      ScamacIdx k = sm->cind[j];
      // current entry is (i, k).
      // try to find corresponding entry (k, i). If sm is symmetric, the "row"ptr is equal to the "column"ptr
      if (ivec[k] >= sm->rptr[k+1]) {// no corresponding element left in this row.
        my_s_pat=false;
        my_s_val=false;
      } else if (sm->cind[ivec[k]] != i) {// by assumption, cind is ordered for each row. Therefore, no corresponding element found.
        my_s_pat=false;
        my_s_val=false;
      } else {// corresponding element found.
        if (sm->valtype == SCAMAC_VAL_REAL) {
          if (sm->val[j] != sm->val[ivec[k]]) {// no value symmetry
            my_s_val=false;
          }
        } else {//SCAMAC_VAL_COMPLEX
          if (sm->val[2*j] != sm->val[2*ivec[k]] || sm->val[2*j+1] != - sm->val[2*ivec[k]+1]) {// conjugate
            my_s_val=false;
          }
        }
        ivec[k]++;
      }
    }
    if (!my_s_pat) {
      break;
    }
  }

  free(ivec);

  if (symm_pattern) {
    *symm_pattern=my_s_pat;
  }
  if (symm_value) {
    *symm_value=my_s_val;
  }

  return SCAMAC_EOK;
}
