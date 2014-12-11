#include "phist_config.h"
#include "phist_typedefs.h"
#include "phist_MultiVector.hpp"
#include "phist_BelosAdapter.hpp"

#ifdef PHIST_HAVE_BELOS
/* this file is just to put some static member variables somewhere */ 
namespace Belos {

#ifdef PHIST_HAVE_SP
int MultiVecTraits<float,::phist::MultiVector<float> >::ierr_;
int MultiVecTraits<s_complex_t, ::phist::MultiVector<s_complex_t> >::ierr_;
#endif
int MultiVecTraits<double, ::phist::MultiVector<double> >::ierr_;
int MultiVecTraits<d_complex_t, ::phist::MultiVector<d_complex_t> >::ierr_;
}
#endif
