//! internal implementation of the orthog function, the regular SUBR(orthog) should be used instead in actual algorithms.

//! which does not randomize the null space (if any). Note that the arguments are renamed and reordered, this should 
//! eventually be adjusted to avoid confusion. Here, W is the orthogonal space, V is orthogonalized against W, and the
//! relation that holds after the subroutine (with Q=V on output) is Q*R1 = V-W*R2, Q'Q=I
void SUBR(orthogrr)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, 
        TYPE(const_sdMat_ptr) WtW_I, TYPE(sdMat_ptr) VtV, _MT_ desiredEps, int maxIter, int* iflag);

//! B-orthogonalization interface, BV is in- and output (contains BQ afterwards)
void SUBR(orthogrrB)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(mvec_ptr) BV,
        TYPE(const_linearOp_ptr) B_op, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, 
        TYPE(const_sdMat_ptr) WtW_I, TYPE(sdMat_ptr) VtV, _MT_ desiredEps, int maxIter, int* iflag);
