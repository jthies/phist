#include "phist_config.h"

#include "phist_precon.h"

namespace phist {

//
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

  static void Create(st::linearOp_t* P, 
        const void* A, ST sigma, const void* B, std::string options, int* iflag)
  {
    NotImplemented(iflag);
    return;
  }
  static void Apply(ST alpha, void const* P, st::mvec_t* X, ST beta, st:mvec_t const* b)
  {
    NotImplemented(iflag);
    return;
  }
  static void ApplyT(ST alpha, void const* P, st::mvec_t* X, ST beta, st:mvec_t const* b)
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
