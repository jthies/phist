#ifndef ESMAC_AUX_H
#define ESMAC_AUX_H

#include "esmac_types.h"

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif
#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX3
#define MAX3(x,y,z) (((x) > MAX(y,z)) ? (x) : MAX(y,z))
#endif
#ifndef MIN3
#define MIN3(x,y,z) (((x) < MIN(y,z)) ? (x) : MIN(y,z))
#endif

/* a multi-index counter */

int * esmac_counter_alloc(int n);
void esmac_counter_reset(int n, int *c);
// increase counter by 1, such that c[i]<=maxc[i]
// returns 0 if counter reached last state
int esmac_counter_step(int n, const int *maxc, int *c);

// return a number larger than n
int esmac_increase_n_somewhat(int n);

// sort (esmac_idx_t) array and eliminates doublets. Returns number of elements after sorting
int esmac_sort_purge_array(int n, esmac_idx_t *arr);

#endif /* ESMAC_AUX_H */
