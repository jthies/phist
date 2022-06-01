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

namespace phist {

//! map owner object
class MapOwner
{

  public:

    //! constructor
    MapOwner(phist_map_ptr v=nullptr){v_=v;}

    //! destructor
    ~MapOwner(){int iflag=0; if (v_!=nullptr) phist_map_delete(v_,&iflag);}

    //! set map pointer
    void set(phist_map_ptr v) {v_=v;}

    //! get map pointer
    phist_map_ptr get() {return v_;}

    //! get map pointer (const)
    phist_const_map_ptr get() const {return v_;}

  private:

    //! wraped map pointer
    phist_map_ptr v_;
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

}//namespace phist

#define PHIST_CLASSFILE_DEF "phist_MemOwner_decl.hpp"
#include "phist_gen_all.h"

#endif
