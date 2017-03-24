/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef KERNELS_PETSC_TYPEDEFS_HPP
#define KERNELS_PETSC_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <petsc.h>
#include <petscmat.h>

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"


template <typename ST>
class Traits
{

public:
  
  //! multi vectors
  struct mvec_t
  {
    Mat v;
    bool is_view = false;
    ST *rawData = NULL;
    phist_const_map_ptr map;
  };

  //! serial dense matrix
  struct sdMat_t
  {
    Mat m;
    bool is_view = false;
    ST* rawData = NULL;
    phist_lidx lda;
    phist_const_comm_ptr comm;
  };

  //! CRS matrices
  struct sparseMat_t
  {
    Mat m;
    phist_const_map_ptr map;
  };

};
  
#endif
