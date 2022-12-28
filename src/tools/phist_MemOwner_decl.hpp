/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag);
extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag);
extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) V, int* iflag);
//! \file phist_MemOwner_decl.hpp
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
//!
//! or, if you do not have an object yet, you may want to use the constructor that actually creates the mvec:
//!
//! mvecOwner<ST> V(map_ptr, ncols, &iflag);
//!
//! and retrieve the pointer using V.get() for calling phist functions.
//!
namespace phist {

//! mvec owner object
template<> class MvecOwner<_ST_>
{

  public:

    //! constructor that takes ownership of an existing pointer if an argument is given
    MvecOwner(TYPE(mvec_ptr) v=nullptr){v_=v;}

    //! constructor that allocates memory
    MvecOwner(phist_const_map_ptr map, int ncols, int *iflag)
    {
      v_=nullptr;
      PHIST_CHK_IERR(SUBR(mvec_create)(&v_, map, ncols, iflag), *iflag);
    }

    //! destructor
    ~MvecOwner(){int iflag=0; if (v_!=nullptr) SUBR(mvec_delete)(v_,&iflag);}

    //! set mvec pointer
    void set(TYPE(mvec_ptr) v) {v_=v;}

    //! get mvec pointer
    TYPE(mvec_ptr) get() {return v_;}

    //! get mvec pointer (const)
    TYPE(const_mvec_ptr) get() const {return v_;}

  private:

    //! wraped mvec pointer
    TYPE(mvec_ptr) v_;

};


//! sdMat owner object
template<> class SdMatOwner<_ST_>
{

  public:

    //! constructor that takes ownership of an existing sdMat_ptr
    SdMatOwner(TYPE(sdMat_ptr) v=nullptr){v_=v;}

    //! constructor that creates a new sdMat and owns it
    SdMatOwner(int nrows, int ncols, phist_const_comm_ptr comm, int* iflag)
    {
      v_=nullptr;
      PHIST_CHK_IERR(SUBR(sdMat_create)(&v_, nrows, ncols, comm, iflag), *iflag);
    }

    //! destructor
    ~SdMatOwner(){int iflag=0; if (v_!=nullptr) SUBR(sdMat_delete)(v_,&iflag);}

    //! set sdMat pointer
    void set(TYPE(sdMat_ptr) v) {v_=v;}

    //! get sdMat pointer
    TYPE(sdMat_ptr) get() {return v_;}

    //! get sdMat pointer (const)
    TYPE(const_sdMat_ptr) get() const {return v_;}

  private:
  
    //! wraped sdMat pointer
    TYPE(sdMat_ptr) v_;
};


//! sparseMat owner object
template<> class SparseMatOwner<_ST_>
{

  public:

    //! constructor
    SparseMatOwner(TYPE(sparseMat_ptr) v=nullptr){v_=v;}

    //! destructor
    ~SparseMatOwner(){int iflag=0; if (v_!=nullptr) SUBR(sparseMat_delete)(v_,&iflag);}

    //! set sparseMat pointer
    void set(TYPE(sparseMat_ptr) v) {v_=v;}

    //! get sparseMat pointer
    TYPE(sparseMat_ptr) get() {return v_;}

    //! get sparseMat pointer (const)
    TYPE(const_sparseMat_ptr) get() const {return v_;}

  private:

    //! wraped sparseMat pointer
    TYPE(sparseMat_ptr) v_;
};

}//namespace phist
