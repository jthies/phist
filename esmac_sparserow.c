#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_aux.h"
#include "esmac_sparserow.h"

esmac_sparserow_real_t * esmac_sparserow_real_alloc() {
  esmac_sparserow_real_t *row = malloc(sizeof *row);
  row->nalloc = 100;
  row->row = malloc(row->nalloc * sizeof *(row->row));
  row->nrow = 0;
  return row;
}

esmac_sparserow_cplx_t * esmac_sparserow_cplx_alloc() {
  esmac_sparserow_cplx_t *row = malloc(sizeof *row);
  row->nalloc = 100;
  row->row = malloc(row->nalloc * sizeof *(row->row));
  row->nrow = 0;
  return row;
}

void esmac_sparserow_real_free(esmac_sparserow_real_t *row) {
  if (row) {
    if (row->row) {free(row->row);}
    free(row);
  }
}

void esmac_sparserow_cplx_free(esmac_sparserow_cplx_t *row) {
  if (row) {
    if (row->row) {free(row->row);}
    free(row);
  }
}

void esmac_sparserow_real_zero(esmac_sparserow_real_t *row) {
  if (row) {
    row->nrow=0;
  }
}

void esmac_sparserow_cplx_zero(esmac_sparserow_cplx_t *row) {
  if (row) {
    row->nrow=0;
  }
}


int esmac_sparserow_real_add(esmac_sparserow_real_t *row, double alpha, esmac_idx_t idx) {
  if (!row) {
    printf("%s: argument row == null. Abort.\n",__func__);
    exit(EXIT_FAILURE);
  }
  if (alpha != 0.0 && idx >=0) {
    (row->nrow)++;
    if (esmac_increase_n_somewhat(row->nrow) > row->nalloc) {
      row->nalloc = esmac_increase_n_somewhat(row->nrow);
      row->row = realloc(row->row, row->nalloc * sizeof *(row->row) );
    }
    row->row[row->nrow-1].val=alpha;
    row->row[row->nrow-1].idx=idx;
  }
  return row->nrow;
}

int esmac_sparserow_cplx_add(esmac_sparserow_cplx_t *row, double complex alpha, esmac_idx_t idx) {
  if (!row) {
    printf("%s: argument row == null. Abort.\n",__func__);
    exit(EXIT_FAILURE);
  }
  if (alpha != 0.0 && idx >=0) {
    (row->nrow)++;
    if (esmac_increase_n_somewhat(row->nrow) > row->nalloc) {
      row->nalloc = esmac_increase_n_somewhat(row->nrow);
      row->row = realloc(row->row, row->nalloc * sizeof *(row->row) );
    }
    row->row[row->nrow-1].val=alpha;
    row->row[row->nrow-1].idx=idx;
  }
  return row->nrow;
}

static int compare_idxreal(const void *a, const void *b) {
  const esmac_idxreal_t * mya = (const esmac_idxreal_t *) a;
  const esmac_idxreal_t * myb = (const esmac_idxreal_t *) b;
  if (mya->idx < myb->idx) {
    return -1;
  } else if (mya->idx > myb->idx) {
    return 1;
  } else {
    return 0;
  }
}


int esmac_sparserow_real_to_idxval(esmac_sparserow_real_t *row, int maxrowlength, int keep_zeros, esmac_idx_t *idx, double *val) {
  if (row->nrow>0) {
    // sort
    qsort(row->row,row->nrow,sizeof *(row->row),compare_idxreal);
    
    // eliminate doublets
    int i,j;    
    i=0; 
    for (j=1;j<row->nrow;j++) {
      if (row->row[i].idx != row->row[j].idx) {
	if (keep_zeros || row->row[i].val != 0.0) {
	  i++;
	}
	if (i<j) {
	  row->row[i] = row->row[j];
	}
      } else {
	row->row[i].val += row->row[j].val;
      }
    } 
    //last one
    if (keep_zeros || row->row[i].val != 0.0) {
      i++;
    }
    row->nrow=i;
    if (idx) {
      if (i>maxrowlength) {
	printf("%s: i>maxrowlength. Abort.\n",__func__);
	exit(EXIT_FAILURE);
      }
      for (i=0;i<row->nrow;i++) {
	idx[i]=row->row[i].idx;
	if (val) {
	  val[i]=row->row[i].val;
	}
      }
    }
    return row->nrow;
  } else {
    return 0;
  }
}

int esmac_sparserow_cplx_to_idxval(esmac_sparserow_cplx_t *row, int maxrowlength, int keep_zeros, esmac_idx_t *idx, double complex *val) {
  if (row->nrow>0) {
    // sort
    qsort(row->row,row->nrow,sizeof *(row->row),compare_idxreal);
    
    // eliminate doublets
    int i,j;    
    i=0; 
    for (j=1;j<row->nrow;j++) {
      if (row->row[i].idx != row->row[j].idx) {
	if (keep_zeros || row->row[i].val != 0.0) {
	  i++;
	}
	if (i<j) {
	  row->row[i] = row->row[j];
	}
      } else {
	row->row[i].val += row->row[j].val;
      }
    } 
    //last one
    if (keep_zeros || row->row[i].val != 0.0) {
      i++;
    }
    row->nrow=i;
    if (idx) {
      if (i>maxrowlength) {
	printf("%s: i>maxrowlength. Abort.\n",__func__);
	exit(EXIT_FAILURE);
      }
      for (i=0;i<row->nrow;i++) {
	idx[i]=row->row[i].idx;
	if (val) {
	  val[i]=row->row[i].val;
	}
      }
    }
    return row->nrow;
  } else {
    return 0;
 }
}

