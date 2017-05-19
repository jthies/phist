#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_aux.h"
#include "esmac_mstate.h"

/* multi-state routines */

esmac_mstate_t * esmac_mstate_alloc(int n) {
  esmac_mstate_t * ms = malloc(sizeof *ms);
  ms->n = n;
  ms->ns = malloc(n * sizeof *(ms->ns));
  ms->ns_total=0;
  ms->qs = malloc(n * sizeof *(ms->qs));
  ms->idx = malloc(n * sizeof *(ms->idx));
  ms->which_dof = malloc(n * sizeof *(ms->which_dof));
  ms->dof = malloc(n * sizeof *(ms->dof));
  int i;
  for (i=0;i<n;i++) {
    ms->ns[i]=0;
    ms->qs[i]=NULL;
    ms->which_dof[i]=ESMAC_DOF_NONE;
    ms->dof[i]=NULL;
  }
  ms->nalloc = 100;
  ms->row = malloc(ms->nalloc * sizeof *(ms->row));
  ms->nrow = 0;
  ms->cnt = malloc(n * sizeof *(ms->cnt));
  ms->maxcnt = malloc(n * sizeof *(ms->maxcnt));
  return ms;
}

void esmac_mstate_free(esmac_mstate_t *ms) {
  if (ms) {
    if (ms->ns) {free(ms->ns);}
    if (ms->qs) {free(ms->qs);}
    if (ms->idx) {free(ms->idx);}
    if (ms->which_dof) {free(ms->which_dof);}
    if (ms->dof) {free(ms->dof);}
    if (ms->row) {free(ms->row);}
    if (ms->cnt) {free(ms->cnt);}
    if (ms->maxcnt) {free(ms->maxcnt);}
    free(ms);
  }
}

// set QS at position 0 <= pos < n, with ns states per QS
int esmac_mstate_set_qs(esmac_mstate_t *ms, int pos, esmac_idx_t ns, esmac_qstate_t *qs) {
  if (0 <= pos && pos < ms->n) {
    ms->qs[pos]=qs;
    ms->ns[pos]=ns;
    // recalculate
    int i;
    ms->ns_total=1;
    for (i=0;i<ms->n;i++) {
      ms->ns_total=ms->ns_total * ms->ns[i];
    }
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}

/*
int esmac_mstate_register_dof_fermions(esmac_mstate_t *ms, int pos, esmac_dof_fermions_t *dof) {
  if (0 <= pos && pos < ms->n) {
    ms->which_dof[pos]=ESMAC_DOF_FERMIONS;
    ms->dof[pos]=dof;
    int info = esmac_mstate_set_qs(ms, pos,esmac_dof_fermions_ns(dof), esmac_dof_fermions_qs(dof) );
    return info;
  } else {
    return ESMAC_EINVAL;
  }
}
*/

int esmac_mstate_register_dof_bosons(esmac_mstate_t *ms, int pos, esmac_dof_bosons_t *dof) {
  if (0 <= pos && pos < ms->n) {
    ms->which_dof[pos]=ESMAC_DOF_BOSONS;
    ms->dof[pos]=dof;
    int info = esmac_mstate_set_qs(ms, pos,esmac_dof_bosons_ns(dof), esmac_dof_bosons_qs(dof) );
    return info;
  } else {
    return ESMAC_EINVAL;
  }
}



/* multistate routines */

// obtain mstate from idx
// implies esmac_mstate_row_zero
int esmac_mstate_decode(esmac_mstate_t *ms, esmac_idx_t idx) {
  if (0<=idx && idx < ms->ns_total) {
    int i;
    for (i=0;i<ms->n;i++) {
      ms->idx[i] = idx % ms->ns[i];
      idx = idx / ms->ns[i];
      if (ms->which_dof[i] != ESMAC_DOF_NONE && ms->dof[i]) {
        switch(ms->which_dof[i]) {
          /*
          case ESMAC_DOF_FERMIONS :
            esmac_dof_fermions_init((esmac_dof_fermions_t *) ms->dof[i], ms->idx[i]);
            break;
            */
          case ESMAC_DOF_BOSONS :
            esmac_dof_bosons_init((esmac_dof_bosons_t *) ms->dof[i], ms->idx[i]);
            break;
          default:
            printf("%s: Unknown DOF\n",__func__);
            exit(EXIT_FAILURE);
            break;
        }
      }
    }
    esmac_mstate_row_zero(ms);
    return 0;
  }  else {
    return ESMAC_EINVAL;
  }
}

// individual idx for qs at pos (with 0 <= idx < 
// error for idx < 0
esmac_idx_t esmac_mstate_qs_idx(esmac_mstate_t *ms, int pos) {
  if (0 <= pos && pos < ms->n) {
    return ms->idx[pos];
  } else {
    return -1;
  }
}

void esmac_mstate_row_zero(esmac_mstate_t *ms) {
  ms->nrow=0;
}

static int compare_idxval(const void *a, const void *b) {
  const esmac_idxval_t *mya = (const esmac_idxval_t *) a;
  const esmac_idxval_t *myb = (const esmac_idxval_t *) b;
  if (mya->idx < myb->idx) {
    return -1;
  } else if (mya->idx > myb->idx) {
    return 1;
  } else {
    return 0;
  }
}

static int compress_row(int n, esmac_idxval_t *row) {
  if (n<=1) {// nothing to do
    return n;
  }
  // sort
  qsort(row,n,sizeof(esmac_idxval_t),compare_idxval);

  // eliminate doubles
  // (keeps ZERO entries)

  int i,j;

  i=0; 
  for (j=1;j<n;j++) {
    if (row[i].idx != row[j].idx) {
      i++;
      if (i<j) {
	row[i] = row[j];
      }
    } else {
      row[i].val = row[i].val + row[j].val;
    }
  }

  return i+1;

}


//add current active qs (as manipulated elsewhere) to output row
// if active part of one qs is empty, nothing is added here
int esmac_mstate_row_add(esmac_mstate_t *ms, double alpha) {
  if (alpha != 0.0) {
    int total=1;
    int i;
    for (i=0;i<ms->n;i++) {
      ms->cnt[i]=0;
      int k;
      k=esmac_qstate_n(ms->qs[i]);
      if (k) {
	ms->maxcnt[i]=k-1;
	total=total*k;
      } else {
	total=0;
	break;
      } 
    }
    if (total) {
      while (1) {
	(ms->nrow)++;
	if (esmac_increase_n_somewhat(ms->nrow) > ms->nalloc) {
	  ms->nalloc = esmac_increase_n_somewhat(ms->nrow);
	  ms->row = realloc(ms->row, ms->nalloc * sizeof *(ms->row) );
	}
	esmac_idx_t accidx=0, idx;
	double accval=1.0, val;
	for (i=(ms->n)-1;i>=0;i--) {
	  esmac_qstate_get(ms->qs[i],ms->cnt[i],&idx,&val); // ignore return value
	  accidx=ms->ns[i]*accidx + idx;
	  accval=accval*val;
	}
	ms->row[ms->nrow-1].idx = accidx;
	ms->row[ms->nrow-1].val = alpha*accval;
	if (!esmac_counter_step(ms->n, ms->maxcnt, ms->cnt)) {
	  break;
	}
      }
    }
    //compress row, eliminate doublets
    ms->nrow = compress_row(ms->nrow,ms->row);
  }
  // reset qstates
  int i;
  for (i=0;i<ms->n;i++) {
    esmac_qstate_reset(ms->qs[i]);
  }
  return 0;
}

int esmac_mstate_row_add_one_qs(esmac_mstate_t *ms, int pos, double alpha) {
  
  if (alpha != 0.0) {
    int total=esmac_qstate_n(ms->qs[pos]);
    if (total) {
      esmac_idx_t accidx = 0, idxscal = 1;
      int i;
      for (i=(ms->n)-1;i>=0;i--) {
        accidx=ms->ns[i]*accidx;
        if (i != pos) {
          accidx = accidx + ms->idx[i];
        }
        if (i < pos) {
          idxscal = idxscal * ms->ns[i];
        }
	    }
      int k;
      for (k=0;k<total;k++) {
        (ms->nrow)++;
        if (esmac_increase_n_somewhat(ms->nrow) > ms->nalloc) {
          ms->nalloc = esmac_increase_n_somewhat(ms->nrow);
          ms->row = realloc(ms->row, ms->nalloc * sizeof *(ms->row) );
        }
        esmac_idx_t idx;
        double val;
        esmac_qstate_get(ms->qs[pos],k,&idx,&val); // ignore return value
        ms->row[ms->nrow-1].idx = accidx + idxscal * idx;
        ms->row[ms->nrow-1].val = alpha*val;
      }
    }
    //compress row, eliminate doublets
    ms->nrow = compress_row(ms->nrow,ms->row);
  }
  // reset qstate
  esmac_qstate_reset(ms->qs[pos]);
  
  return 0;
}

/* get one row of a sparse matrix from the multi-state */
int esmac_mstate_row_to_idxval(esmac_mstate_t *ms, int maxrowlength, int keep_zeros, esmac_idx_t *idx, double *val) {
  if (idx && val) {
    if (keep_zeros && ms->nrow > maxrowlength) {
      return ESMAC_ESHORTROW;
    }
    if ( ms->nrow > maxrowlength) { // count precisely, it may just fit
      int i, n=0;
      for (i=0;i<ms->nrow;i++) {
	if (ms->row[i].val != 0.0) {
	  n++;
	}
      }
      if (n >  maxrowlength) {
	return ESMAC_ESHORTROW;
      }
    }
    int i, n=0;
    for (i=0;i<ms->nrow;i++) {
      if (keep_zeros || ms->row[i].val != 0.0) {
	idx[i]=ms->row[n].idx;
	val[i]=ms->row[n].val;
	n++;
      }
    }
    return n;
  } else {
    if (keep_zeros) {
      return ms->nrow;
    } else { // count precisely
      int i, n=0;
      for (i=0;i<ms->nrow;i++) {
	if (ms->row[i].val != 0.0) {
	  n++;
	}
      }
      return n;
    }
  }
}

