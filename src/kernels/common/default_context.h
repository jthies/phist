/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_INTERNAL_DEFAULT_CONTEXT_HPP
#define PHIST_INTERNAL_DEFAULT_CONTEXT_HPP

#include "phist_config.h"
#include "phist_void_aliases.h"
#include "phist_kernels.h"

namespace phist 
{
  namespace internal 
  {

    class default_context 
    {
      public:
      
      default_context(phist_const_map_ptr row_map=NULL,
                      phist_const_map_ptr col_map=NULL,
                      phist_const_map_ptr range_map=NULL,
                      phist_const_map_ptr domain_map=NULL);

      ~default_context();

      phist_const_map_ptr row_map;
      phist_const_map_ptr col_map;
      phist_const_map_ptr range_map;
      phist_const_map_ptr domain_map;
    };
    
    typedef std::map<const void*,default_context*> contextMap;

    static contextMap contextCollection;
    
    //! delete any context objects associated with the memory location A
    void delete_default_context(const void* A);
    //! retrieve a context associated with the memory location A, or create a
    //! new one and associate it with A.
    default_context* get_default_context(const void* A);
  
  }
}


#endif
