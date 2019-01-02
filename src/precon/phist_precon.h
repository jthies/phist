/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_precon.h 
//! \brief Preconditioning interface

#ifndef PHIST_PRECON_H
#define PHIST_PRECON_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"
#include "phist_kernels.h"



#include "phist_enums.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup precon precon: Preconditioning interface
//!@{

//! \brief data structure used internally for storing preconditioner data
//!
//! this will be the A_ field in the linearOp_t, but the
//! user should only work with linearOp instead of this one directly
typedef struct {
  //! \brief identifies the preconditioner type (e.g. IFPACK)
  //!
  //! This is used to internally call the correct create/delete/apply etc functions
  phist_Eprecon type_;
  void const* A_; //!< pointer to matrix A
  void const* B_; //!< pointer to matrix B
  void const *Vkern_; //!< pointer to vector space approximating the kernel of A-sigma*B 
  void const *BVkern_; //!< pointer to vector space approximating the kernel of A-sigma*B 
  void* P_; //!< pointer to preconditioning object
} phist_internal_precon;

#define PHIST_CLASSFILE_DEF "phist_precon_decl.h"
#include "phist_gen_all.h"

//!@}

#ifdef __cplusplus
}
#endif


#endif
