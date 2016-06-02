#ifndef PHIST_GHOST_PRIVATE_H
#define PHIST_GHOST_PRIVATE_H

#ifndef PHIST_HAVE_MPI
typedef int MPI_Comm;
const int MPI_COMM_WORLD=0;
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


//! map object for ghost kernel library
//!                                                                                        
//! in the petra object model that this abstract kernel interface follows, a map describes 
//! the distribution of data among processors. It is used to e.g. define how vector entries
//! or matrix rows are distributed among MPI processes and required to create a vector.    
//!                                                                                        
//! In ghost this is handled differently, the 'context' is the main object for describing  
//! communication patterns, and it is owned by a sparse matrix. In order to create a vector
//! we need the vtraits_t object, which also knows about number of vectors and             
//! data type. So we define our own struct with a pointer to the context object and a temp-
//! late for cloning the vtraits in mvec_create.
//!
//! A map can be obtained in two ways in phist: by map_create or by asking a matrix for its
//! row/range/col/domain map. In the former case we will return a uniformly distributed map
//! with a newly created context. Vectors (mvecs) based on such an object can not be multi-
//! plied by a matrix. The maps obtained from a sparseMat are presently identical (row/range/
//! col/domain). This is because ghost keeps the communication buffer for spmv in the vector,
//! not the matrix. So any object that can be multiplied by a matrix requires the "column map",
//! and it is not possible to include a vector permutation or redistribution in the spmv, as
//! is allowed by the Petra object model.
//!
//! The only way to convert an mvec from one map into another is by calling the function
//! mvec_to_mvec. 
typedef struct ghost_map
  {
  ghost_context* ctx;
  ghost_densemat_traits vtraits_template;
  } ghost_map;


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

    //! this calls ghost_context_create with the given arguments and retries with a proc_weight of 1 if there are empty partitions
    //! If the matrix source is not NULL (i.e. the matrix is read from disk or created from a row function) 
    //! gnrows=gncols=0 should be given.
    void context_create(ghost_context **context, ghost_gidx gnrows, ghost_gidx gncols, 
        ghost_context_flags_t flags, void *matrixSource, ghost_sparsemat_src srcType, ghost_mpi_comm comm, double proc_weight, int* iflag);
    
    //! private helper function to create a vtraits object
    ghost_densemat_traits phist_default_vtraits();
    
    //! A Garbage collection for maps as they need to be recreated dynamically with GHOST.
    
    //! Each time someone requests a map from an mvec or sparseMat, a new object is created
    //! and addedto an internal 'collection'. Since maps do not contain data but only pointers
    //! we currently never delete these objects until the end of the run.
    //!
    //! \todo we could use a hash map instead and associate each map with the context it links to.
    //! Reference-counting could then be used to see if a context is no longer needed and delete it.
    //! 
    //! Example: sparseMat is created       => context created @X
    //!          sparseMat_get_domain_map   => map[X] created with ref count 1
    //!          mvec_create(&V,map[X],...) => ref count of map[X] incremented
    //!          sparseMat_delete           => ref count of map[X] decremented
    //!          mvec_delete(V)             => ref count of map[X] decremented to 0, map and context deleted
    //! See issue #177
    class MapGarbageCollector
    {

      public:
      
        //!
        ghost_map* new_map(const void* p);
        //!
        void delete_maps(void* p);

      private:
        //!
        typedef std::map<const void*, std::vector<ghost_map*> > MapCollection;
        //!
        MapCollection maps_;
    };
  
    static MapGarbageCollector mapGarbageCollector;
  
  }
}

#endif
