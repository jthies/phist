#include "phist_config.h"
#include "phist_typedefs.h"
#include "phist_MultiVector.hpp"
#include "phist_BelosAdapter.hpp"

#ifdef PHIST_HAVE_BELOS
/* this file is just to put some static member variables somewhere */ 
namespace Belos {

#ifdef PHIST_HAVE_SP
int MultiVecTraits<float,::phist::MultiVector<float> >::iflag_;
int MultiVecTraits<s_complex_t, ::phist::MultiVector<s_complex_t> >::iflag_;
#endif
int MultiVecTraits<double, ::phist::MultiVector<double> >::iflag_;
int MultiVecTraits<d_complex_t, ::phist::MultiVector<d_complex_t> >::iflag_;
}

namespace phist 
{
#ifdef PHIST_HAVE_SP
  int MultiVector<float>::countObjects=0;
  int MultiVector<s_complex_t>::countObjects=0;
#endif
  int MultiVector<double>::countObjects=0;
  int MultiVector<d_complex_t>::countObjects=0;
}
#endif
