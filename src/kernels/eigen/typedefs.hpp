/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef KERNELS_EIGEN_TYPEDEFS_HPP
#define KERNELS_EIGEN_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "phist_typedefs.h"
#include "phist_ScalarTraits.hpp"


template <typename ST>
class Traits
{

public:
  
  //! multi vectors
  struct mvec_t
  {
#ifdef PHIST_MVECS_ROW_MAJOR
    typedef Eigen::Matrix<ST, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Eigen_MVec;
#else
    typedef Eigen::Matrix<ST, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Eigen_MVec;
#endif
    Eigen_MVec v_storage;
    typedef Eigen::Ref<Eigen_MVec> Eigen_Ref;
    Eigen_Ref v = v_storage;
    phist_const_map_ptr map;
  };

  //! serial dense matrix
  struct sdMat_t
  {
    typedef Eigen::Matrix<ST, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Eigen_SdMat;
    Eigen_SdMat m_storage;
    typedef Eigen::Ref<Eigen_SdMat> Eigen_Ref;
    Eigen_Ref m = m_storage;
    phist_const_comm_ptr comm;
  };

  //! CRS matrices
  struct sparseMat_t
  {
    typedef Eigen::SparseMatrix<ST> Eigen_SparseMat;
    Eigen_SparseMat m;
    phist_const_map_ptr map;
  };

};
  
#endif
