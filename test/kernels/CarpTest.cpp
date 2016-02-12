#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#include "kernels/phist_kernels.h"
#include "kernels/phist_kernel_flags.h"
#include "KernelTestWithSparseMat.h"
#include "KernelTestWithVectors.h"

#include <complex>

using namespace ::testing;

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

#define MATNAME MATNAME_BENCH3D_8_A1
#define _BASENAME_ CarpTest
/* the tests are only executed if the carp kernel is implemented,
if the setup routine returns -99, further tests will not be run
*/
#warning "Reenable again for real cases"
#define DISABLE_TESTCASE

#define _N_ 512
#define _M_ 1
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 4
#include "../phist_typed_test_gen.h"

#define _N_ 512
#define _M_ 7
#include "../phist_typed_test_gen.h"
