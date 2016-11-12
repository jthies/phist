//! This file contains simple classes to wrap mvec, sdMat and sparseMat
//! pointers in C++ objects which will take care of their deletion when
//! destroyed themselves. There is no reference counting mechanism, but
//! the wrapper objects could be put into std::shared_ptr's to achieve 
//! automatic memory management.
//!
//! A typical use of these classes is:
//! 
//! {
//! TYPE(mvec_ptr) tmp_V;
//! PHIST_CHK_IERR(SUBR(mvec_create)(&tmp_V,map,nvec,iflag),*iflag);
//! mvecOwner<ST> ownV(tmpV);
//! PHIST_CHK_IERR(... some other calls ...)
//! ...
//! wherever the program leaves the scope, tmp_V is deleted without
//! the user having to call SUBR(mvec_delete)
//! }

//!
template<> class MvecOwner<_ST_>
{

  public:

    //!
    MvecOwner(TYPE(mvec_ptr) v=NULL){v_=v;}

    //!
    ~MvecOwner(){int iflag=0; if (v_!=NULL) SUBR(mvec_delete)(v_,&iflag);}
    
    //!
    void set(TYPE(mvec_ptr) v) {v_=v;}

  private:
  
    //!
    TYPE(mvec_ptr) v_;
    
};


//!
template<> class SdMatOwner<_ST_>
{

  public:
  
    //!
    SdMatOwner(TYPE(sdMat_ptr) v=NULL){v_=v;}

    //!
    ~SdMatOwner(){int iflag=0; if (v_!=NULL) SUBR(sdMat_delete)(v_,&iflag);}

    //!
    void set(TYPE(sdMat_ptr) v) {v_=v;}

  private:
  
    //!
    TYPE(sdMat_ptr) v_;
};


//!
template<> class SparseMatOwner<_ST_>
{

  public:

    //!
    SparseMatOwner(TYPE(sparseMat_ptr) v=NULL){v_=v;}

    //!
    ~SparseMatOwner(){int iflag=0; if (v_!=NULL) SUBR(sparseMat_delete)(v_,&iflag);}

    //!
    void set(TYPE(sparseMat_ptr) v) {v_=v;}

  private:
  
    //!
    TYPE(sparseMat_ptr) v_;
};

