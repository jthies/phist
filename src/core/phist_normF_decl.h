
//! \defgroup core Core functionality used by algorithms
//@{

//!\name Frobenius norm for sdMats
void SUBR(sdMat_normF)(TYPE(const_sdMat_ptr) M, _MT_ *f, int* iflag);

//!\name Frobenius norm for mvecs
void SUBR(mvec_normF)(TYPE(const_mvec_ptr) V, _MT_ *f, int* iflag);

//@}


