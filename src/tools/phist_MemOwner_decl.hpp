/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
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

//! mvec owner object
template<> class MvecOwner<_ST_>
{

  public:

    //! constructor
    MvecOwner(TYPE(mvec_ptr) v=NULL){v_=v;}

    //! destructor
    ~MvecOwner(){int iflag=0; if (v_!=NULL) SUBR(mvec_delete)(v_,&iflag);}
    
    //! set mvec pointer
    void set(TYPE(mvec_ptr) v) {v_=v;}

  private:
  
    //! wraped mvec pointer
    TYPE(mvec_ptr) v_;
    
};


//! sdMat owner object
template<> class SdMatOwner<_ST_>
{

  public:
  
    //! constructor
    SdMatOwner(TYPE(sdMat_ptr) v=NULL){v_=v;}

    //! destructor
    ~SdMatOwner(){int iflag=0; if (v_!=NULL) SUBR(sdMat_delete)(v_,&iflag);}

    //! set sdMat pointer
    void set(TYPE(sdMat_ptr) v) {v_=v;}

  private:
  
    //! wraped sdMat pointer
    TYPE(sdMat_ptr) v_;
};


//! sparseMat owner object
template<> class SparseMatOwner<_ST_>
{

  public:

    //! constructor
    SparseMatOwner(TYPE(sparseMat_ptr) v=NULL){v_=v;}

    //! destructor
    ~SparseMatOwner(){int iflag=0; if (v_!=NULL) SUBR(sparseMat_delete)(v_,&iflag);}

    //! set sparseMat pointer
    void set(TYPE(sparseMat_ptr) v) {v_=v;}

  private:
  
    //! wraped sparseMat pointer
    TYPE(sparseMat_ptr) v_;
};

