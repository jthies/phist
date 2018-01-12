/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "./default_context.h"

namespace phist 
{
  namespace internal 
  {

    default_context::default_context(phist_const_map_ptr row_map_in,
                                             phist_const_map_ptr col_map_in,
                                             phist_const_map_ptr range_map_in,
                                             phist_const_map_ptr domain_map_in)
    :
      row_map{row_map_in},
      col_map{col_map_in},
      range_map{range_map_in},
      domain_map{domain_map_in}
    {
    }

    default_context::~default_context()
    {
      // do nothing - we don't own the maps!
    }

    void delete_default_context(const void* A)
    {
      contextMap::iterator it=contextCollection.find(A);
      if (it!=contextCollection.end())
      {
        delete it->second;
        contextCollection.erase(A);
      }
    }
    
    default_context* get_default_context(const void* A)
    {
      default_context* ctx=NULL;
      contextMap::iterator it=contextCollection.find(A);
      if (it!=contextCollection.end())
      {
        ctx=it->second;
      }
      else
      {
        ctx=new default_context();
        contextCollection[A]=ctx;
      }
      return ctx;
    }
  }
}

// create a context object defining the shape of a matrix by three maps.
// The row map determines the distribution of rows of the matrix over the processes.                       
// The range and domain map determine the size and distribution of y and x, resp. in y=A*x, and therefore the
// shape of the sparse matrices created with this context. Note that it is typically not necessary to create a
// context a priori, the more common case is creating a sparseMat, obtaining its context and using it to create
// another ("compatible") sparseMat. This gives the kernel library the freedom to apply permutations and load-
// balancing when creating a new matrix.
extern "C" void phist_context_create(phist_context_ptr* vctx, 
                          phist_const_map_ptr row_map, 
                          phist_const_map_ptr col_map,
                          phist_const_map_ptr range_map,
                          phist_const_map_ptr domain_map, 
                          int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  *vctx=(phist_context_ptr)(new phist::internal::default_context
        (row_map,col_map, range_map,domain_map));
}

// delete a context created by context_create
extern "C" void phist_context_delete(phist_context_ptr vctx, int* iflag)
{
  phist::internal::default_context* ctx=(phist::internal::default_context*)vctx;
  delete ctx;
}
