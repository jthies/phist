#if PHIST_OUTLEV>=PHIST_TRACE
#define TRACE_PHISTMV_MEM
#endif

namespace phist {

//! a C++ wrapper for ghost_vec_t so that  
//! we can add a destructor. If we don't   
//! do this, the Trilinos RCP concept won't
//! work and memory leaks will occur.      
//! In order to allow creating an RCP like 
//! this:                                  
//!                                        
//! Teuchos::RCP<ghost_vec_t> v_ptr        
//!      =phist::rcp(v);                   
//!                                        
//! we derive MultiVector from ghost_vec_t. Be-
//! ware that the 'member functions' of    
//! ghost_vec_t are NULL, however, unless  
//! you pass the object through            
//! ghost_vecCreate(). The clean way of    
//! using this wrapper is to always pass   
//! in a pointer to a complete ghost_vec_t 
//! (after ghost_createVec) and use the    
//! get() method to access it.             
//!                                        
//! The ownMem flag is similar to the flag 
//! you pass to an rcp, if true the vector 
//! is destroyed along with the wrapper,   
//! otherwise just the wrapper is deleted. 
//!                                        
//! JT 20.11.2013: disabling the feature   
//! that MultiVector is derived from ghost_vec 
//! because it leads ot confusion and bugs.
//! Use phist::rcp(ghost_vec_t*) to get an 
//! RCP<MultiVector> instead, and phist::ref2ptr
//! in order to get a ghost_vec_t* from a  
//! MultiVector object (or whatever the underly-
//! ing kernel lib uses as vector format).
template<>
class MultiVector< _ST_ >
{
#include "phist_std_typedefs.hpp"

  public:

  //!
  MultiVector(mvec_ptr_t v_in, bool ownMem)
    {
    v_=v_in;
    ownMem_=ownMem;
    
#ifdef TRACE_PHISTMV_MEM
    myID=countObjects++;
    PHIST_OUT(PHIST_INFO,"### Create MultiVector #%d, ownMem=%d\n",myID,ownMem);
#endif
    }
  
  //!
  virtual ~MultiVector()
    {
#ifdef TRACE_PHISTMV_MEM
    if (ownMem_)
      {
      PHIST_OUT(PHIST_INFO,"### Delete MultiVector #%d\n",myID);
      }
    else
      {
      PHIST_OUT(PHIST_INFO,"### Delete view MultiVector #%d\n",myID);
      }
#endif
    if (ownMem_)
      {
      PHIST_TCHK_IERR(SUBR(mvec_delete)(v_,&iflag_),iflag_);
      this->v_=NULL;
      }
    }

  //!
  mvec_ptr_t get()
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for mvec_t";
    }

  //!
  const_mvec_ptr_t get() const
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for mvec_t";
    }
  
protected:

  //! disallow default constructor
  MultiVector()
    {
    v_=NULL;
    }
  
  //! disallow copy constructor
  MultiVector(const MultiVector& v)
    {
    (void)v;//unused
    v_=NULL;
    }
  
  //! disallow assignment operator
  MultiVector& operator=(const MultiVector& v)
    {
    MultiVector *w = new MultiVector(v);
    return *w;
    }
  
  
  //! the wrapped object
  mvec_ptr_t v_;

  //! are we allowed to delete the vector?
  bool ownMem_;

  // give each newly created object a label so we can track where they are destroyed 
  static int countObjects;
  
  // internally used error flag
  int iflag_;

  // label of this object
  int myID;
};

} //namespace phist

