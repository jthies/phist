#ifndef ESMAC_LUT_H
#define ESMAC_LUT_H

#include "esmac_types.h"

/*   _L_ook_U_p _T_ables  (LUT)
 * 
 */

esmac_idx_t esmac_lut_construct(int ineq, int n, int m, int s, esmac_idx_t **lut);
esmac_idx_t esmac_lut_encode(int ineq, int n, int s, esmac_idx_t *cnt, const int *x);
     int esmac_lut_decode(int ineq, int n, int s, esmac_idx_t *cnt, esmac_idx_t idx, int *x);
     
/* for fermions or spins (sz=const): fixed s after construction, ineq = 0, m =1 (i.e, x[i] = 0,1) */

esmac_idx_t esmac_lut_onezero_construct(int n, int s, esmac_idx_t **lut);
esmac_idx_t esmac_lut_onezero_encode(int n, int s, esmac_idx_t *cnt, const int *x);
     int esmac_lut_onezero_decode(int n, int s, esmac_idx_t *cnt, esmac_idx_t idx, int *x);

#endif /* ESMAC_LUT_H */
