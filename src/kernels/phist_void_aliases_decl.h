/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \addtogroup void_alias
//!@{
    
//! opaque multi-vector object
typedef struct {} TYPE(mvec);

//! opaque pointer to multi-vector objects
typedef TYPE(mvec)* TYPE(mvec_ptr);

//! opaque pointer to const multi-vector objects
typedef const TYPE(mvec)* TYPE(const_mvec_ptr);

//! opaque small dense matrix object
typedef struct {} TYPE(sdMat);

//! opaque pointer to small dense matrix objects
typedef TYPE(sdMat)* TYPE(sdMat_ptr);

//! opaque pointer to const small dense matrix objects
typedef const TYPE(sdMat)* TYPE(const_sdMat_ptr);

//! opaque CRS matrix
typedef struct {} TYPE(sparseMat);

//! opaque pointer to CRS matrix objects
typedef TYPE(sparseMat)* TYPE(sparseMat_ptr);

//! opaque pointer to const CRS matrix objects
typedef const TYPE(sparseMat)* TYPE(const_sparseMat_ptr);

//!@}