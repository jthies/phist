//! apply a transformation matrix M to a given search space V and AV, BV and the projection H=V'AV
void SUBR(transform_searchSpace)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) M, bool generalizedEigenproblem, int* iflag);

//! apply a transformation matrices M to given search spaces V and W, BV 
//! and the projections H=W'V, H_A=W'AV:
//!
//! V   <- V*MV
//! W   <- W*MW
//! H   <- MW'*H*MV
//! H_A <- MW'*H*MV
void SUBR(transform_searchSpaceHarmonic)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) W, TYPE(mvec_ptr) BV,
        TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) H_A, TYPE(sdMat_ptr) MV, TYPE(sdMat_ptr) MW,
        bool generalizedEigenproblem, int* iflag);
