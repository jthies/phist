#ifndef PHIST_PRECON_TRAITS_HPP
#define PHIST_PRECON_TRAITS_HPP

#include "phist_config.h"

namespace phist {

//! traits class for making preconditioners usable in phist programs. The user
//! should not use this class directly but the C wrapper defined in phist_precon.h

//! default implementation of the traits class that allows us to interface
//! with our own and third-party preconditioners. The default implementation
//! just gives error messages about missing template specialization.
template<typename ST,precon_t PT>
class PreconTraits
{
  typedef ScalarTraits<ST> st;

  static void NotImplemented(int* iflag)
  {
    PHIST_SOUT(PHIST_ERROR("class PreconTraits is missing a specialization for data type %s and preconditioner type %s\n",
                        st::type_char(), precon2str(PT));
    *iflag=PHIST_NOT_IMPLEMENTED;
  }

  static void Create(void** P, 
        const void* A, ST sigma, const void* B, const char* options, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void Delete(void* P, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void Apply(ST alpha, void const* P, st::mvec_t const* X, ST beta, st:mvec_t* Y)
  {
    NotImplemented(iflag);
    return;
  }
  static void ApplyT(ST alpha, void const* P, st::mvec_t const* X, ST beta, st:mvec_t* Y)
  {
    NotImplemented(iflag);
    return;
  }
  static void ApplyShifted((ST alpha, const void* P, ST const * sigma,
          st::mvec_t const* X, ST beta,  st::mvec_t* Y, int* iflag);
  {
    NotImplemented(iflag);
    return;
  }
};

}

// depending on the installation, include source files with template specializations
# if defined(PHIST_KERNEL_LIB_EPETRA)
# include "tpl/phist_Ifpack_def.hpp"
# include "tpl/phist_ML_def.hpp"
# include "tpl/phist_MUELU_def.hpp"
# elif defined(PHIST_KERNEL_LIB_TPETRA
# include "tpl/phist_Ifpack2_def.hpp"
# include "tpl/phist_Amesos2_def.hpp"
# include "tpl/phist_MUELU_def.hpp"
# endif

#endif
