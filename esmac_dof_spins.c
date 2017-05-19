#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "esmac_constants.h"
#include "esmac_lut.h"
#include "esmac_dof_spins.h"

esmac_dof_spins_t * esmac_dof_spins_alloc(int n_sites) {
  esmac_dof_spins_t * mydof = malloc(sizeof *mydof);

  mydof->n_sites=n_sites;
  mydof->n_up=-1; // no constraint

  //mydof->ns = esmac_lut_construct(0, n_sites, 1, n_spins, &(mydof->cnt));
  mydof->ns = 1U << n_sites;
  mydof->cnt=NULL;
  
  return mydof;
}

esmac_dof_spins_t * esmac_dof_spinssz_alloc(int n_sites, int n_up) {
  esmac_dof_spins_t * mydof = malloc(sizeof *mydof);

  mydof->n_sites=n_sites;
  mydof->n_up=n_up; // no constraint

  //mydof->ns = esmac_lut_construct(0, n_sites, 1, n_spins, &(mydof->cnt));
  mydof->ns = esmac_lut_onezero_construct(n_sites, n_up, &(mydof->cnt));
  
  return mydof;
}



void esmac_dof_spins_free(esmac_dof_spins_t * dof) {
  if (dof) {
    if (dof->cnt) {free(dof->cnt);}
    free(dof);
  }
}

esmac_idx_t esmac_dof_spins_ns(const esmac_dof_spins_t * dof) {
  return dof->ns;
}

esmac_rep_spins_t * esmac_rep_spins_alloc(const esmac_dof_spins_t *dof) {
  esmac_rep_spins_t * myrep = malloc(dof->n_sites * sizeof *myrep);
  return myrep;
}

void esmac_rep_spins_free(esmac_rep_spins_t * rep) {
  if (rep) {free(rep);}
}

void esmac_rep_spins_copy(const esmac_dof_spins_t *dof, const esmac_rep_spins_t * rep, esmac_rep_spins_t * repcpy) {
  if (rep && repcpy) {
    memcpy(repcpy, rep, dof->n_sites * sizeof *rep );
  }
}

int esmac_spins_decode(const esmac_dof_spins_t * dof, esmac_idx_t idx, esmac_rep_spins_t *rep) {
  if (0 <= idx && idx < dof->ns) {
   // esmac_lut_decode(0, dof->n_sites, dof->n_spins, dof->cnt, idx, dof->rep);
    if (dof->n_up<0) {
      int i;
      for (i=0;i<dof->n_sites;i++) {
        if (idx & 1) {// bitwise
          rep[dof->n_sites-i-1]=1;
        } else {
          rep[dof->n_sites-i-1]=0;
        }
        idx = idx >> 1;
      }
    } else {
      esmac_lut_onezero_decode(dof->n_sites, dof->n_up, dof->cnt, idx, rep);
    }
    return 0;
  } else {
    return ESMAC_EINVAL;
  }
}

esmac_idx_t esmac_spins_encode(const esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep) {
  if (dof->n_up<0) {
    esmac_idx_t idx=0;
    int i;
    for (i=0;i<dof->n_sites;i++) {
      idx = idx << 1;
      if (rep[i]) {
        idx++;
      }
    }
    return idx;
  } else {
    return esmac_lut_onezero_encode(dof->n_sites, dof->n_up, dof->cnt, rep);
  }
}

/* operators */

double esmac_op_spins_Sz   (esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep, int i) {
  if (0 <= i && i < dof->n_sites) {
    if (rep[i]) {
      return 0.5;
    } else {
      return -0.5;
    }
  } else {
    return 0.0;
  }
}

double esmac_op_spins_SzSz (esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    if (rep[i] == rep[j]) {
      return 0.25;
    } else {
      return -0.25;
    }
  } else {
    return 0.0;
  }
}

double esmac_op_spins_SzSum(esmac_dof_spins_t * dof, const esmac_rep_spins_t *rep) {
  if (dof->n_up >=0) {
    return ((double) dof->n_up - 0.5*dof->n_sites);
  } else {
    int i;
    double val = 0.0;
    for (i=0;i<dof->n_sites;i++) {
      if (rep[i]) {
        val = val + 0.5;
      } else {
        val = val - 0.5;
      }
    }
    return val;
  }
}

double esmac_op_spins_Sp   (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i) {
  assert( dof->n_up <0 );
  if (0 <= i && i < dof->n_sites) {
    if (rep[i]) {
      return 0.0;
    } else {
      rep[i]=1;
      return 1.0;
    }
  } else {
    return 0.0;
  }
} 

double esmac_op_spins_Sm   (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i) {
  assert( dof->n_up <0 );
  if (0 <= i && i < dof->n_sites) {
    if (rep[i]) {
      rep[i]=0;
      return 1.0;
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
} 

double esmac_op_spins_Sx   (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i) {
  assert( dof->n_up <0 );
  if (0 <= i && i < dof->n_sites) {
    rep[i]=1-rep[i];
    return 0.5;
  } else {
    return 0.0;
  }
}


double esmac_op_spins_SpSp (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j) {
  assert( dof->n_up <0 );
  if (i == j) {return 0.0;}
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    if (rep[i] || rep[j]) {
      return 0.0;
    } else {
      rep[i]=1;
      rep[j]=1;
      return 1.0;
    }
  } else {
    return 0.0;
  }
}

double esmac_op_spins_SpSm (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j) {
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    if (i==j) {
      if (rep[i]) {
        return 1.0;
      } else {
        return 0.0;
      }
    } else {
      if (!rep[i] && rep[j]) {
        rep[i]=1;
        rep[j]=0;
        return 1.0;
      } else {
        return 0.0;
      }
    }
  } else {
    return 0.0;
  }
}

double esmac_op_spins_SmSm (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j) {
  assert( dof->n_up <0 );
  if (i == j) {return 0.0;}
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    if (rep[i] && rep[j]) {
      rep[i]=0;
      rep[j]=0;
      return 0.0;
    } else {
      return 1.0;
    }
  } else {
    return 0.0;
  }
}

double esmac_op_spins_SxSx (esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j) {
  assert( dof->n_up <0 );
  if (i == j) {return 0.25;}
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    rep[i]=1-rep[i];
    rep[j]=1-rep[j];
    return 0.25;
  } else {
    return 0.0;
  }
}

double esmac_op_spins_Sflip(esmac_dof_spins_t * dof,       esmac_rep_spins_t *rep, int i, int j) {
  if (i == j) {return 0.0;}
  if (0 <= i && i < dof->n_sites && 0 <= j && j < dof->n_sites) {
    if (rep[i] != rep[j]) {
      rep[i]=1-rep[i];
      rep[j]=1-rep[j];
      return 0.5;
    } else {
      return 0.0;
    }
  } else {
    return 0.0;
  }
}
