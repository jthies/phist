#ifndef ESMAC_CONSTANTS_H
#define ESMAC_CONSTANTS_H

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


// error codes
#define ESMAC_EINVAL 1
#define ESMAC_ESHORTROW -7
#define ESMAC_ERANGE -132
#define ESMAC_ENULL -323

// length of names in ESMAC
#define ESMAC_NAME_LENGTH 128

// identifiers for degrees of freedom
#define ESMAC_DOF_NONE 0
#define ESMAC_DOF_FERMIONS 1
#define ESMAC_DOF_BOSONS 2

// flags for esmac_matrix_info_t
#define ESMAC_VAL_REAL 1
#define ESMAC_VAL_COMPLEX 2
#define ESMAC_NONE 1
#define ESMAC_SYMMETRIC 2
#define ESMAC_HERMITIAN 3

// boundary conditions
#define ESMAC_OBC 1
#define ESMAC_PBC 2

// parameter types
#define ESMAC_PAR_INT 1
#define ESMAC_PAR_DOUBLE 2
#define ESMAC_PAR_BOOL 3


#endif /* ESMAC_CONSTANTS_H */
