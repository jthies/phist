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

    //! this calls ghost_context_create with the given arguments and retries with a proc_weight of 1 if there are empty partitions
    //! If the matrix source is not NULL (i.e. the matrix is read from disk or created from a row function) 
    //! gnrows=gncols=0 should be given.
    void context_create(ghost_context **context, ghost_gidx gnrows, ghost_gidx gncols, 
        ghost_context_flags_t flags, void *matrixSource, ghost_sparsemat_src srcType, ghost_mpi_comm comm, double proc_weight, int* iflag);
    
    //! private helper function to create a vtraits object
    ghost_densemat_traits phist_default_vtraits();
    

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
class ghost_map
{

  public:

  ghost_context* ctx;
  ghost_densemat_permutation *perm_local;
  ghost_densemat_permutation *perm_global;
  // blueprints for how to create a densemat (mvec) object from this map
  ghost_densemat_traits vtraits_template;
  // blueprints for how to create a sparsemat object from this map
  ghost_sparsemat_traits mtraits_template;
  
  //! create a new map which gets a context, a permutation type (defining which permutation object, row/col/none of the
  //! context should be applied to an mvec based on the map, and two flags:
  //! own_ctx=true means that when the map is deleted, the context should be deleted, too.
  //! if own_ctx=true and own_perm=false, the permutation objects of the context are *not* deleted,
  //! this is used to share permutation objects among multiple contexts/maps.
  ghost_map(ghost_context* ctx_in, ghost_densemat_permuted pt=NONE, bool own_ctx=false, bool own_perm=true)
  {
    perm_local=NULL;
    perm_global=NULL;
    ctx=ctx_in;
    ownContext_=own_ctx;
    ownPermutations_=own_perm;

    vtraits_template=phist_default_vtraits();
    mtraits_template=(ghost_sparsemat_traits)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
    // we set C=sigma=-1 to indicate that the map has not been used for creating any matrix.
    // The functions that actually generate the maps (mvec/sparseMat_get*map) will overrule 
    // this setting, for instance, if you create a sparseMat from a map obtained by an mvec,
    // the sparseMat can not be sorted globally or locally anymore because it has to be
    mtraits_template.C=-2;
    mtraits_template.sortScope=-2;
    mtraits_template.flags=GHOST_SPARSEMAT_DEFAULT;
    if (pt==COLUMN)
    {
      // we need to set a few things by hand to make sure we can
      // check compatibility of maps
      int me;
      ghost_rank(&me,ctx->mpicomm);
      vtraits_template.gnrows=ctx->gnrows;
      vtraits_template.nrows=ctx->lnrows[me];
    }
    else
    {
      vtraits_template.gnrows=ctx->gncols;
      // NOTE: if we have non-square matrices, the number of
      //       local rows of the vector must be determined  
      //       differently, for now we leave it at its default
      //       value (0 or -1, I guess) and let ghost decide.
      //       such a map can currently not be comapred with
      //       maps_compatible, however!
      if (ctx->gncols==ctx->gnrows)
      {
        int me;
        ghost_rank(&me,ctx->mpicomm);
        vtraits_template.nrows=ctx->lnrows[me];
      }
    }
        
    if (pt!=NONE && (ctx->perm_local || ctx->perm_global))
    {
      vtraits_template.permutemethod=pt;
      vtraits_template.flags|=GHOST_DENSEMAT_PERMUTED;
  
      perm_local=new ghost_densemat_permutation;
      perm_global=new ghost_densemat_permutation;
      perm_local->perm=NULL;
      perm_local->invPerm=NULL;
      perm_global->perm=NULL;
      perm_global->invPerm=NULL;
    }
    if (pt==ROW)
    {
      if (ctx->perm_local) 
      {
        perm_local->perm=ctx->perm_local->perm;
        perm_local->invPerm=ctx->perm_local->invPerm;
      }
      if (ctx->perm_global) 
      {
        perm_global->perm=ctx->perm_global->perm;
        perm_global->invPerm=ctx->perm_global->invPerm;
      }
    }
    else if (pt==COLUMN)
    {
      if (ctx->perm_local) 
      {
        perm_local->perm=ctx->perm_local->colPerm;
        perm_local->invPerm=ctx->perm_local->colInvPerm;
      }
      if (ctx->perm_global) 
      {
        perm_global->perm=ctx->perm_global->colPerm;
        perm_global->invPerm=ctx->perm_global->colInvPerm;
      }
    }
  }
  
  ~ghost_map()
  {
    if (ownContext_)
    {
      if (!ownPermutations_)
      {
        // someone else has to take care of this,
        // we're sharing the data with another context
        ctx->perm_local=NULL;
        ctx->perm_global=NULL;
      }
      ghost_context_destroy(ctx);
    }
    // note: the context owns data structures, these
    // perm objects only have pointers to data in the context,
    // so we can always allocate the objects without too much 
    // overhead
    if (perm_local) delete perm_local;
    if (perm_global) delete perm_global;
  }
  
  private:
  
  //! if true, the context will be destroyed when the map is (i.e. in phist_map_delete)
  bool ownContext_;

  //! if true, the context's perm objects will be destroyed when the context is in phist_map_delete
  bool ownPermutations_;
  
};


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
      
        //! create a new map and associate it with the object pointed to by p. If there are already maps
        //! associated with pointer p, return the first one found *if reuse_if_exists==true*.
        ghost_map* new_map(const void* p, ghost_context* ctx=NULL, ghost_densemat_permuted pt=NONE, 
                bool own_ctx=false, bool own_perm=true, bool reuse_if_exists=false);
        //! associate an existing map with the object pointed to by p
        void add_map(const void* p, ghost_map* m);
        //! delete all maps associated with the object pointed to by p
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
