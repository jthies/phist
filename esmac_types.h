#ifndef ESMAC_TYPES_H
#define ESMAC_TYPES_H

#include <stdint.h>
#include <inttypes.h>

#include "esmac_config.h"


#ifdef ESMAC_INDEX_TYPE_int64 
  typedef int64_t esmac_idx_t;
  #define ESMACPRIDX PRId64
#elif defined ESMAC_INDEX_TYPE_int32
  typedef int32_t esmac_idx_t;
  #define ESMACPRIDX PRId32
#else
  typedef int esmac_idx_t;
  #error "Unknown ESMAC_INDEX_TYPE"
#endif

//#define esmac_idx_t long int
//#define ESMACPRIDX "ld"

typedef struct {
  esmac_idx_t nrow;
  esmac_idx_t ncol;
  int maxrowlength;
  int valtype; // ESMAC_VAL_REAL or ESMAC_VAL_COMPLEX
  int symmetry; // ESMAC_NONE or ESMAC_SYMMETRIC or ESMAC_HERMITIAN
} esmac_matrix_info_t;


#endif /* ESMAC_TYPES_H */
