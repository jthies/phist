#ifndef PHIST_DRIVER_UTILS_H
#define PHIST_DRIVER_UTILS_H

#ifndef _ST_
#error "this file should be included after a phist_gen_X header"
#endif

typedef _ST_ ST;
typedef _MT_ MT;
typedef TYPE(mvec_ptr) mvec_ptr_t;
typedef TYPE(const_mvec_ptr) const_mvec_ptr_t;

typedef TYPE(sdMat_ptr) sdMat_ptr_t;
typedef TYPE(const_sdMat_ptr) const_sdMat_ptr_t;

typedef TYPE(crsMat_ptr) crsMat_ptr_t;
typedef TYPE(const_crsMat_ptr) const_crsMat_ptr_t;

typedef TYPE(op_ptr) op_ptr_t;
typedef TYPE(const_op_ptr) const_op_ptr_t;

//! read matrix from some supported file format

//! auto-detects the file type by looking at the file extension
//! and calls the appropriate kernel routine to read the matrix
void SUBR(crsMat_read)(TYPE(crsMat_ptr)* A, char* filename, int* ierr);

#endif
