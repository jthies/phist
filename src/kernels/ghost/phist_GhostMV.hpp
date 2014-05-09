#ifndef BELOS_GHOST_MV_HPP
#define BELOS_GHOST_MV_HPP

#include "phist_config.h"
#include "phist_macros.h"
#include "ghost.h"

#if PHIST_OUTLEV>=PHIST_TRACE
#define TRACE_GHOSTMV_MEM
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
//! JT 20.11.2013: disabling the feature   
//! that GhostMV is derived from ghost_vec 
//! because it leads ot confusion and bugs.
//! Use phist::rcp(ghost_vec_t*) to get an 
//! RCP<GhostMV> instead, and phist::ref2ptr
//! in order to get a ghost_vec_t* from a  
//! GhostMV object (or whatever the underly-
//! ing kernel lib uses as vector format).
class GhostMV
  {
  public:
  //!
  GhostMV(ghost_densemat_t* v_in, bool ownMem)
    {
    v_=v_in;
    ownMem_=ownMem;
    
#ifdef TRACE_GHOSTMV_MEM
    myID=countObjects++;
    PHIST_OUT(PHIST_INFO,"### Create GhostMV #%d, ownMem=%d",myID,ownMem);
#endif
    }
  
  //!
  virtual ~GhostMV()
    {
#ifdef TRACE_GHOSTMV_MEM
    if (ownMem_)
      {
      PHIST_OUT(PHIST_INFO,"### Delete GhostMV #%d",myID);
      }
    else
      {
      PHIST_OUT(PHIST_INFO,"### Delete view GhostMV #%d",myID);
      }
#endif
    if (ownMem_)
      {
      this->get()->destroy(this->get());
      this->v_=NULL;
      }
    }

  //!
  ghost_densemat_t* get()
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_densemat_t";
    }

  //!
  const ghost_densemat_t* get() const
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_densemat_t";
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
  ghost_densemat_t* v_;

  //! are we allowed to delete the vector?
  bool ownMem_;

  // give each newly created object a label so we can track where they are destroyed 
  static int countObjects;

#if 1
//#if PHIST_OUTLEV>=PHIST_TRACE  
  // label of this object
  int myID;
#endif
  };

} //namespace phist

#endif
