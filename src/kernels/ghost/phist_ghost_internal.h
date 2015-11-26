#ifndef PHIST_GHOST_PRIVATE_H
#define PHIST_GHOST_PRIVATE_H

#ifndef PHIST_HAVE_MPI
typedef int MPI_Comm;
const int MPI_COMM_WORLD=0;
#endif


#ifndef __cplusplus
#error "this is a C++ header"
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
//! late for cloning the vtraits in mvec_create
typedef struct ghost_map_t
  {
  ghost_context_t* ctx;
  ghost_densemat_traits_t vtraits_template;
  ghost_permutation_t *permutation;
  } ghost_map_t;


    //! small helper function to preclude integer overflows (ghost allows 64 bit local indices, 
    //! but we don't right now)
    template<typename idx_t>
    int check_local_size(idx_t& i)
    {
      if (i>std::numeric_limits<lidx_t>::max())
      {
        return -1;
      }
      return 0;
    }


    //! set reasonable default parameters for SELL-C-sigma sparse matrix format in GHOST
    void get_C_sigma(int* C, int* sigma, int flags, MPI_Comm comm);

    //! private helper function to create a vtraits object
    ghost_densemat_traits_t phist_default_vtraits();
    
    //! A Garbage collection for maps as they need to be recreated dynamically with GHOST
    class MapGarbageCollector
    {
      public:
      
        //!
        ghost_map_t* new_map(const void* p);
        //!
        void delete_maps(void* p);

      private:
        //!
        typedef std::map<const void*, std::vector<ghost_map_t*> > MapCollection;
        //!
        MapCollection maps_;
    };
  
    static MapGarbageCollector mapGarbageCollector;
  
  }
}

#endif
