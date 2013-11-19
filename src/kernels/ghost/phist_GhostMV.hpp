#ifndef BELOS_GHOST_MV_HPP
#define BELOS_GHOST_MV_HPP

#include "phist_macros.h"
#include "ghost.h"

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
//! we derive GhostMV from ghost_vec_t. Be-
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
class GhostMV: public ghost_vec_t
  {
  public:
  //!
  GhostMV(ghost_vec_t* v_in, bool ownMem)
    {
    v_=v_in;
    ownMem_=ownMem;
    
#if PHIST_OUTLEV>=PHIST_TRACE
    PHIST_OUT(PHIST_TRACE,"+++ Create GhostMV #%d, ownMem=%d",myID,ownMem);
    myID=countObjects++;
#endif
    }
  
  //!
  virtual ~GhostMV()
    {
    if (ownMem_)
      {
      this->get()->destroy(this->get());
      }
#if PHIST_OUTLEV>=PHIST_TRACE
    if (ownMem_)
      {
      PHIST_OUT(PHIST_TRACE,"+++ Delete GhostMV #%d",myID);
      }
    else
      {
      PHIST_OUT(PHIST_TRACE,"+++ Delete view GhostMV #%d",myID);
      }
#endif
    }

  //!
  ghost_vec_t* get()
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_vec_t";
    }

  //!
  const ghost_vec_t* get() const
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_vec_t";
    }
  
protected:

  //! disallow default constructor
  GhostMV()
    {
    v_=NULL;
    }
  
  //! disallow copy constructor
  GhostMV(const GhostMV& v)
    {
    (void)v;//unused
    v_=NULL;
    }
  
  //! disallow assignment operator
  GhostMV& operator=(const GhostMV& v)
    {
    GhostMV *w = new GhostMV(v);
    return *w;
    }
  
  
  //! the wrapped object
  ghost_vec_t* v_;

  //! are we allowed to delete the vector?
  bool ownMem_;

#if PHIST_OUTLEV>=PHIST_TRACE  
  // give each newly created object a label so we can track where they are destroyed 
  static int countObjects;
  // label of this object
  int myID;
#endif  
  };

} //namespace phist

#endif
