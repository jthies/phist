/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_LUT_H
#define SCAMAC_LUT_H

#include <stdbool.h>

#include "scamac_include.h"

/*   _L_ook_U_p _T_ables  (LUT)
 *
 */

// return total number of possibilites.
// return negative value, if overflow occurs.
ScamacErrorCode scamac_lut_construct(bool ineq, int n, int m, int s, ScamacIdx * ns, ScamacIdx ** lut);

ScamacIdx scamac_lut_encode(bool ineq, int n, int s, ScamacIdx * cnt, const int * x);
void scamac_lut_decode(bool ineq, int n, int s, ScamacIdx * cnt, ScamacIdx idx, int * x);

/* for fermions or spins (sz=const): fixed s after construction, ineq = 0, m =1 (i.e, x[i] = 0,1) */

ScamacErrorCode scamac_lut_onezero_construct(int n, int s, ScamacIdx * ns, ScamacIdx ** lut);
ScamacIdx scamac_lut_onezero_encode(int n, int s, ScamacIdx * cnt, const int *x);
void scamac_lut_onezero_decode(int n, int s, ScamacIdx * cnt, ScamacIdx idx, int * x);

#endif /* SCAMAC_LUT_H */
