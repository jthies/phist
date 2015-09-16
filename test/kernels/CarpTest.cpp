#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"
//#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "kernels/phist_kernel_flags.h"
#include "KernelTestWithVectors.h"
#include "../tools/MatrixIO.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <complex>

using namespace ::testing;

extern "C" {

int phist_Sidfunc(ghost_gidx_t row, ghost_lidx_t* len, ghost_gidx_t* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((float*)vals)[0]=1.0f;
  return 0;
}
int phist_Didfunc(ghost_gidx_t row, ghost_lidx_t* len, ghost_gidx_t* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((double*)vals)[0]=1.0;
  return 0;
}
int phist_Cidfunc(ghost_gidx_t row, ghost_lidx_t* len, ghost_gidx_t* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((std::complex<float>*)vals)[0]=(std::complex<float>)1.0f;
  return 0;
}
int phist_Zidfunc(ghost_gidx_t row, ghost_lidx_t* len, ghost_gidx_t* cols, void* vals, void *arg)
{
  *len=1;
  cols[0]=row;
  ((std::complex<double>*)vals)[0]=(std::complex<double>)1.0;
  return 0;
}

}//extern "C"

#define _N_ 512
#define _MATNAME_ "BENCH3D-8-A1"
#define _NV_ 1

/* the tests are only executed if the carp kernel is implemented,
if the setup routine returns -99, further tests will not be run
*/
#ifdef PHIST_HAVE_SP
# define CLASSNAME SCarpTest_512_1
# include "phist_gen_s.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME DISABLED_CCarpTest_512_1
# include "phist_gen_c.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME
#endif

#define CLASSNAME DCarpTest_512_1
#include "phist_gen_d.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME DISABLED_ZCarpTest_512_1
#include "phist_gen_z.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#undef _NV_
#define _NV_ 4

#ifdef PHIST_HAVE_SP

# define CLASSNAME SCarpTest_512_4
# include "phist_gen_s.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME DISABLED_CCarpTest_512_4
# include "phist_gen_c.h"
# include "CarpTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DCarpTest_512_4
#include "phist_gen_d.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME DISABLED_ZCarpTest_512_4
#include "phist_gen_z.h"
#include "CarpTest_def.hpp"

#undef CLASSNAME
#undef _NV_
#define _NV_ 7

#ifdef PHIST_HAVE_SP

# define CLASSNAME SCarpTest_512_7
# include "phist_gen_s.h"
# include "CarpTest_def.hpp"
# undef CLASSNAME

# define CLASSNAME DISABLED_CCarpTest_512_7
# include "phist_gen_c.h"
# include "CarpTest_def.hpp"
#undef CLASSNAME

#endif

#define CLASSNAME DCarpTest_512_7
#include "phist_gen_d.h"
#include "CarpTest_def.hpp"
#undef CLASSNAME

#define CLASSNAME DISABLED_ZCarpTest_512_7
#include "phist_gen_z.h"
#include "CarpTest_def.hpp"
