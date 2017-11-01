/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_AUX_H
#define SCAMAC_AUX_H

#include "scamac_include.h"

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

int * scamac_counter_alloc(int n);
void scamac_counter_reset(int n, int *c);
// increase counter by 1, such that c[i]<=maxc[i]
// returns 0 if counter reached last state
int scamac_counter_step(int n, const int *maxc, int *c);

// return a number larger than n
int scamac_increase_n_somewhat(int n);

// sort (ScamacIdx) array and eliminates doublets. Returns number of elements after sorting
int scamac_sort_purge_array(int n, ScamacIdx *arr);

#endif /* SCAMAC_AUX_H */
