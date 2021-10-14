/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_GHOST_PRIVATE_H
#define PHIST_GHOST_PRIVATE_H

#ifndef PHIST_HAVE_MPI
typedef int MPI_Comm;
// note: ghost #defines MPI_COMM_WORLD, I think...
#ifndef MPI_COMM_WORLD
const int MPI_COMM_WORLD=0;
#endif
#endif


#ifndef __cplusplus
#error "this is a C++ header"
#endif

#ifndef DOXYGEN
#include <map>
#endif

namespace phist 
{
  //! some helper functions and objects that should not be exposed to phist users
  namespace ghost_internal 
  {

    //! this calls ghost_context_create with the given arguments and retries with a proc_weight of 1 if there are empty partitions
    //! If the matrix source is not NULL (i.e. the matrix is read from disk or created from a row function) 
    //! gnrows=gncols=0 should be given.
    void context_create(ghost_context **context, ghost_gidx gnrows, ghost_gidx gncols, 
        ghost_context_flags_t flags, void *matrixSource, ghost_sparsemat_src srcType, ghost_mpi_comm comm, double proc_weight, int* iflag);
    
    //! private helper function to create a vtraits object
    ghost_densemat_traits default_vtraits();
    
    //! small helper function to preclude integer overflows (ghost allows 64 bit local indices, 
    //! but we don't right now)
    template<typename idx_t>
    int check_local_size(idx_t& i)
    {
      if (i>std::numeric_limits<phist_lidx>::max())
      {
        return -1;
      }
      return 0;
    }
    
    // get the most promising repartitioning flag (Zoltan, Scotch or nothing)
    // depending on what is available in the ghost installation. if PHIST_OUTLEV>=outlev,
    // a statement will be printed to indicate the selection.
    int get_perm_flag(int iflag, int outlev);

    //! set reasonable default parameters for SELL-C-sigma sparse matrix format in GHOST
    void get_C_sigma(int* C, int* sigma, int flags, MPI_Comm comm);

    //! returns the local partition size based on a benchmark, the benchmark is run once when
    //! this function is called first, subsequently the same weight will be returned unless
    //! you set force_value>0. If force_vale<=0 is given, the benchmark is run anyway and the
    //! newly measured value is returned in this and subsequent calls.
    //! The resulting double can be passed as 'weight' parameter
    //! to ghost_context_create (used in phist_map_create and sparseMat construction routines)
    double get_proc_weight(double force_value=-1.0);

  }
}

#endif
