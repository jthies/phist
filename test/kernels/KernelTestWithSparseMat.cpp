#include "KernelTestWithSparseMat.h"

using namespace phist::testing;

extern "C" {

int phist_Sidfunc(ghost_gidx row, ghost_lidx* len, ghost_gidx* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((float*)vals)[0]=1.0f;
  return 0;
}
int phist_Didfunc(ghost_gidx row, ghost_lidx* len, ghost_gidx* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((double*)vals)[0]=1.0;
  return 0;
}
int phist_Cidfunc(ghost_gidx row, ghost_lidx* len, ghost_gidx* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((std::complex<float>*)vals)[0]=(std::complex<float>)1.0f;
  return 0;
}
int phist_Zidfunc(ghost_gidx row, ghost_lidx* len, ghost_gidx* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((std::complex<double>*)vals)[0]=(std::complex<double>)1.0;
  return 0;
}

}//extern "C"

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
