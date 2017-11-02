// this is what the C++ interface should roughly look like
#include "phist_config.h"
#include "phist_kernels.h"
#include "phist_ScalarTraits.hpp"

namespace phist {

class PhistError: public std::exception 
{
};


class PhistWarning: public std::exception 
{
};

template<typename ST>
class KernelTraits
{
  public: 
  
    typedef typename ST::MissingImplementationOfScalarTraitsClass error;

};

// specializations

template<>
class KernelTraits<double>
{

  typedef double ST;
  typedef ScalarTraits<ST> st;
  typedef ScalarTraits<ST>::magn_t MT;
  typedef TypeTraits<ST> tt;

  // original function comment
  void mvec_create(tt::mvec_ptr* V, const_map_ptr map, int* iflag)
  {
    phist_Dmvec_create(V,map,iflag);
    if (*iflag<0) throw Error(phist_retcode2str(*iflag));
    else if (*iflag>0) throw Warning(phist_retcode2str(*iflag));
  }
};

} // namespace phist
