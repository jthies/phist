/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_SAFEINT_H
#define SCAMAC_SAFEINT_H

#include "scamac_include.h"

// return a+b (or a*b) for a,b >=0 if no overflow occurs,
// otherwise return negative value (signals overflow)
ScamacIdx scamac_safe_add (ScamacIdx a, ScamacIdx b);
ScamacIdx scamac_safe_mult(ScamacIdx a, ScamacIdx b);

/*
// c = binomial(a \over\ b)
ScamacErrorCode scamac_safe_binomial(ScamacIdx a, ScamacIdx b, ScamacIdx *c);
*/

#endif /* SCAMAC_SAFEINT_H */
