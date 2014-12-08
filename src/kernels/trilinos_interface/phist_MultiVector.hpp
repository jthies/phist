#ifndef PHIST_MULTIVECTOR_HPP
#define PHIST_MULTIVECTOR_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"

#include "phist_ScalarTraits.hpp"

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
template<typename __ST>
class MultiVector
{
// all instances of this class that are actually used
// will be specialized below
};

} //namespace phist

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_MultiVector_def.hpp"
#include "phist_gen_c.h"
#include "phist_MultiVector_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_MultiVector_def.hpp"
#include "phist_gen_z.h"
#include "phist_MultiVector_def.hpp"
#endif
