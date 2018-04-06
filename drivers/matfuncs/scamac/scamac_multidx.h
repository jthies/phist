/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_MULTIDX_H
#define SCAMAC_MULTIDX_H

#include "scamac_include.h"

typedef struct {
  int n;
  /* i-th index runs from 0 to ni[i]-1 */
  ScamacIdx *ni;
  /* product of ni[0]*ni[1]* ... */
  ScamacIdx *niprod;
  /* number of differents index entries */
  ScamacIdx nidx;
} scamac_multidx_st;

ScamacErrorCode scamac_multidx_alloc(int n, scamac_multidx_st ** midx);
ScamacErrorCode scamac_multidx_free(scamac_multidx_st *midx);
ScamacErrorCode scamac_multidx_set(scamac_multidx_st *midx, int pos, ScamacIdx ni);

ScamacIdx scamac_multidx_nidx(const scamac_multidx_st *midx);

ScamacIdx scamac_multidx_decode(const scamac_multidx_st *midx, int pos, ScamacIdx idx);

ScamacIdx scamac_multidx_upd(const scamac_multidx_st *midx, int pos, ScamacIdx idx_at_pos, ScamacIdx idx);

#endif /* SCAMAC_MULTIDX_H */
