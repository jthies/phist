#ifndef KERNELS_MAGMA_TYPEDEFS_HPP
#define KERNELS_MAGMA_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "magma.h"

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"


template <typename ST>
class Traits
{

public:
  
  //! multi vectors
  typedef void mvec_t;

  //! serial dense matrix
  typedef void sdMat_t;


  //! CRS matrices
  typedef void sparseMat_t;

};
  
#endif
