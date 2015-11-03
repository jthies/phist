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
  struct mvec_t
  {
    lidx_t n, nvec, stride;
    bool is_view;
    ST* gpuData;
    ST* cpuData;
    const_map_ptr_t map;
  };

  //! serial dense matrix
  struct sdMat_t
  {
    lidx_t nrows, ncols, stride;
    bool is_view;
    ST* gpuData;
    ST* cpuData;
  };

  //! CRS matrices
  typedef void sparseMat_t;

};
  
#endif
