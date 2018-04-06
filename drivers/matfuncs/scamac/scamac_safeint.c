#include "scamac_safeint.h"

ScamacIdx scamac_safe_add (ScamacIdx a, ScamacIdx b) {
  if ( (a < 0) || (b < 0) ) {
    return -1;
  } else if ( a > (SCAMACIDXMAX - b) ) {
    return -2;
  } else {
    return a+b;
  }
}

ScamacIdx scamac_safe_mult(ScamacIdx a, ScamacIdx b) {
  if ( (a < 0) || (b < 0) ) {
    return -1;
  } else if ( b == 0) {
    return 0;
  } else if ( a > (SCAMACIDXMAX / b) ) {
    return -2;
  } else {
    return a*b;
  }
}

/*
ScamacErrorCode scamac_safe_binomial(ScamacIdx a, ScamacIdx b, ScamacIdx *c) {
  if ( (a < 0) || (b < 0) ) {
    return SCAMAC_ERANGE;
  } else {
    reScamacIdx
  }
}
*/

