#include "phist_macros.h"

//! \addtogroup kernels
//@{

/*! @file phist_kernels_common_decl.h 
    @brief simple derived functions implemented for all kernel libs centrally

*/

#ifdef __cplusplus
extern "C" {
#endif


//! \name getting data from objects
//@{

//! get global sparseMat size (number of rows) \ingroup crsmat
void SUBR(sparseMat_global_nrows)(TYPE(sparseMat_ptr) A, gidx_t* s, int* iflag);

//! get global sparseMat size (number of columns) \ingroup crsmat
void SUBR(sparseMat_global_ncols)(TYPE(sparseMat_ptr) A, gidx_t* s, int* iflag);

//! retrieve local length of the vectors in V \ingroup mvec
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* iflag);

//! retrieve global length of the vectors in V \ingroup mvec
void SUBR(mvec_global_length)(TYPE(const_mvec_ptr) V, gidx_t* len, int* iflag);

//! retrieve the comm used for MPI communication in V \ingroup mvec
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* iflag);

//@}

//! y[i]=alpha*(A*x+shift*x) + beta*y
void SUBR(sparseMat_times_mvec_add_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ shift, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);

#ifdef __cplusplus
} // extern "C"
#endif
