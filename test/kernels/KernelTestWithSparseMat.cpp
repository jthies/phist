#include "KernelTestWithSparseMat.h"

using namespace phist::testing;

// available matrices
const char* MatNameEnumToStr(MATNAME_ENUM i)
{
  switch(i)
  {
    case MATNAME_speye:
      return "speye";
    case MATNAME_spzero:
      return "spzero";
    case MATNAME_sprandn:
      return "sprandn";
    case MATNAME_sprandn_nodiag:
      return "sprandn_nodiag";
    case MATNAME_spshift:
      return "spshift";
    case MATNAME_jadaTestMat:
      return "jadaTestMat";
    case MATNAME_symmMat:
      return "symmMat";
    case MATNAME_sprandsym:
      return "sprandsym";
    case MATNAME_IDFUNC:
      return "idfunc";
    case MATNAME_BENCH3D_8_A1:
      return "BENCH3D-8-A1";
    default:
      return "unknown matrix index";
  }
}
