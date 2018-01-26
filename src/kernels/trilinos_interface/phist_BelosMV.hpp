/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_BELOS_MV_HPP
#define PHIST_BELOS_MV_HPP

#include "phist_config.h"
#include "phist_macros.h"
#include "phist_types.hpp"
#include "phist_kernels.hpp"

#ifdef PHIST_HAVE_TEUCHOS
#include "Teuchos_RCP.hpp"
#endif

#if PHIST_OUTLEV>=PHIST_TRACE
#define TRACE_BELOSMV_MEM
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
//! we derive BelosMV from ghost_vec_t. Be-
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
//! that BelosMV is derived from ghost_vec 
//! because it leads ot confusion and bugs.
//! Use phist::rcp(ghost_vec_t*) to get an 
//! RCP<BelosMV> instead, and phist::ref2ptr
//! in order to get a ghost_vec_t* from a  
//! BelosMV object (or whatever the underly-
//! ing kernel lib uses as vector format).
template<typename Scalar>
class BelosMV
{
  public:
  
  typedef typename phist::types<Scalar>::mvec_ptr mvec_ptr;
  typedef typename phist::types<Scalar>::const_mvec_ptr const_mvec_ptr;
  typedef phist::kernels<Scalar> kt;
  
  //!
  BelosMV(mvec_ptr v_in, bool ownMem)
  {
    v_=v_in;
    ownMem_=ownMem;
    
#ifdef TRACE_BELOSMV_MEM
    myID=countObjects++;
    PHIST_OUT(PHIST_INFO,"### Create BelosMV #%d, ownMem=%d\n",myID,ownMem);
#else
    myID=-1;
#endif
  }
  
  //!
  virtual ~BelosMV()
  {
#ifdef TRACE_BELOSMV_MEM
    if (ownMem_)
    {
      PHIST_OUT(PHIST_INFO,"### Delete BelosMV #%d\n",myID);
    }
    else
    {
      PHIST_OUT(PHIST_INFO,"### Delete view BelosMV #%d\n",myID);
    }
#endif
    if (ownMem_)
    {
      int iflag=0;
      kt::mvec_delete(this->get(),&iflag);
      this->v_=nullptr;
    }
  }

  //!
  mvec_ptr get()
    {
    if (v_!=nullptr) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_densemat";
    }

  //!
  const_mvec_ptr get() const
    {
    if (v_!=nullptr) 
      {
      return v_;
      }
    throw "invalid wrapper object for ghost_densemat";
    }
  
protected:

  //! disallow default constructor
  BelosMV()=delete;
  
  //! disallow copy constructor
  BelosMV(const BelosMV<Scalar>& v)=delete;
  
  //! disallow assignment operator
  BelosMV<Scalar>& operator=(const BelosMV<Scalar>& v)=delete;
  
  
  //! the wrapped object
  mvec_ptr v_;

  //! are we allowed to delete the vector?
  bool ownMem_;

  // give each newly created object a label so we can track where they are destroyed 
  static int countObjects;

  // label of this object for TRACE_BELOSMV_MEM
  int myID;
};

#ifdef PHIST_HAVE_TEUCHOS

//! this function is used to easily create the wrapper object from a raw pointer,
// similar to Teuchos::rcp. You do have to explicitly give the scalar type, though,
// because all phist objects are passed around as void pointers. So phist::mvec_rcp<double>(V,true/false)
// instead of just rcp(V,true/false)
template<typename Scalar>
Teuchos::RCP<BelosMV<Scalar> > mvec_rcp(typename phist::types<Scalar>::mvec_ptr V, bool own_mem)
{
  return Teuchos::rcp(new BelosMV<Scalar>(V,own_mem),true);
}

//!
template<typename Scalar>
Teuchos::RCP<const BelosMV<Scalar> > mvec_rcp(typename phist::types<Scalar>::const_mvec_ptr V, bool own_mem)
{
  return Teuchos::rcp(new BelosMV<Scalar>(const_cast<typename phist::types<Scalar>::mvec_ptr>(V),false),true);
}

#endif /* PHIST_HAVE_TEUCHOS */

template<typename ST> int BelosMV<ST>::countObjects=0;

} //namespace phist

#endif
