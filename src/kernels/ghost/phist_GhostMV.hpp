/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef BELOS_GHOST_MV_HPP
#define BELOS_GHOST_MV_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"
#include <ghost.h>

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
  GhostMV(ghost_densemat* v_in, bool ownMem)
    {
    v_=v_in;
    ownMem_=ownMem;
    
#ifdef TRACE_GHOSTMV_MEM
    myID=countObjects++;
    PHIST_OUT(PHIST_INFO,"### Create GhostMV #%d, ownMem=%d\n",myID,ownMem);
#else
    myID=-1;
#endif
    }
  
  //!
  virtual ~GhostMV()
    {
#ifdef TRACE_GHOSTMV_MEM
    if (ownMem_)
      {
      PHIST_OUT(PHIST_INFO,"### Delete GhostMV #%d\n",myID);
      }
    else
      {
      PHIST_OUT(PHIST_INFO,"### Delete view GhostMV #%d\n",myID);
      }
#endif
    if (ownMem_)
      {
      ghost_densemat_destroy(this->get());
      this->v_=NULL;
      }
    }

  //!
  ghost_densemat* get()
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_densemat";
    }

  //!
  const ghost_densemat* get() const
    {
    if (v_!=NULL) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_densemat";
    }
  
protected:

  //! disallow default constructor
  GhostMV();
  
  //! disallow copy constructor
  GhostMV(const GhostMV& v);
  
  //! disallow assignment operator
  GhostMV& operator=(const GhostMV& v);
  
  
  //! the wrapped object
  ghost_densemat* v_;

  //! are we allowed to delete the vector?
  bool ownMem_;

  // give each newly created object a label so we can track where they are destroyed 
  static int countObjects;

  // label of this object for TRACE_GHOSTMV_MEM
  int myID;
  };

} //namespace phist

#endif
