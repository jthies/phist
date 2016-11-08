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

      phist_const_map_ptr range_map;
      phist_const_map_ptr domain_map;
      phist_const_map_ptr row_map;
    };

    static std::map<const void*,default_context*> contextCollection;
    
    void delete_default_context(const void* A)
    {
      if (contextCollection[A]!=NULL)
      {
        delete contextCollection[A];
        contextCollection.erase(A);
      }
    }
  
  }
}


#endif
