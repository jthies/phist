#ifdef __cplusplus
extern "C" {
#endif

typedef struct TYPE(carpData)
{
  TYPE(mvec_ptr) diagA_;
  TYPE(mvec_ptr) rowScaling_;
  TYPE(mvec_ptr) xLoc_; // for importing the X vector
} TYPE(carpData);

//! setup CARP data structures for a matrix A. The data structures
//! are independent of the shift s in an operator A-sB.
void SUBR(carp_create)(TYPE(const_crsMat_ptr) A, TYPE(carpData)** carpData, int* ierr);

//! delete carp data structure
void SUBR(carp_delete)(TYPE(carpData)* carpData, int* ierr);

//! forward/backward sweep of Kaczmarz/CARP algorithm (SSOR sweep on the normal equations),
//! with matrix A-sigma[j]B applied to vector column j (if B==NULL we assume B=I).         
void SUBR(carp_fb)(TYPE(carpData)* carpData, TYPE(const_crsMat_ptr) A, 
        TYPE(const_mvec_ptr) B, _ST_ const * sigma, 
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) X, int* ierr);

#ifdef __cplusplus
} //extern "C"
#endif


