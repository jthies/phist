/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_MEM_OWNER_HPP
#define PHIST_MEM_OWNER_HPP

//! \file phist_MemOwner.hpp
//!
//! This file contains simple classes to wrap mvec, sdMat and sparseMat
//! pointers in C++ objects which will take care of their deletion when
//! destroyed themselves. There is no reference counting mechanism, but
//! the wrapper objects could be put into std::shared_ptr's to achieve 
//! automatic memory management.
//!
//! A typical use of these classes is:
//! 
//!     {
//!     TYPE(mvec_ptr) tmp_V;
//!     PHIST_CHK_IERR(SUBR(mvec_create)(&tmp_V,map,nvec,iflag),*iflag);
//!     mvecOwner<ST> ownV(tmpV);
//!     PHIST_CHK_IERR(... some other calls ...)
//!     ...
//!     wherever the program leaves the scope, tmp_V is deleted without
//!     the user having to call SUBR(mvec_delete)
//!     }

#ifndef DOXYGEN
#include "phist_kernels.h"
#endif

//! map owner object
class MapOwner
{

  public:

    //! constructor
    MapOwner(phist_map_ptr m=nullptr){m_=m;}

    //! destructor
    ~MapOwner(){int iflag=0; if (m_!=nullptr) phist_map_delete)(m_,&iflag);}

    //! set map pointer
    void set(map_ptr m) {m_=m;}

  private:

    //! wraped map pointer
    phist_map_ptr m_;
};

//! mvec owner object
template<typename T> class MvecOwner
{
  typedef typename T::MissingImplementationOfMvecOwnerClass error;
};

//! sdMat owner object
template<typename T> class SdMatOwner
{
  typedef typename T::MissingImplementationOfSdMatOwnerClass error;
};

//! sparseMat owner object
template<typename T> class SparseMatOwner
{
  typedef typename T::MissingImplementationOfSparseMatOwnerClass error;
};

#define PHIST_CLASSFILE_DEF "phist_MemOwner_decl.hpp"
#include "phist_gen_all.h"

#endif
