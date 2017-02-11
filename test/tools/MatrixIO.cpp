#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"

#include "MatrixIO.h"
#include <iostream>

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_ScalarTraits.hpp"

namespace phist
{
  namespace testing
  {

    int getSparseMatCreateFlag(int N, int NV)
    {
      int flag=0;
      if (N>=100)
      {
        if (NV>1)
        {
          flag=PHIST_SPARSEMAT_OPT_BLOCKSPMVM;
        }
        else
        {
          flag=PHIST_SPARSEMAT_OPT_SINGLESPMVM;
        }
      }
      if (flag) flag|=PHIST_SPARSEMAT_PERM_LOCAL;
#if PHIST_OUTLEV<4
      flag|=PHIST_SPARSEMAT_QUIET;
#endif
      return flag;
    }
  }
}
#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "MatrixIO_def.hpp"

#include "phist_gen_c.h"
#include "MatrixIO_def.hpp"

#endif

#include "phist_gen_d.h"
#include "MatrixIO_def.hpp"

#include "phist_gen_z.h"
#include "MatrixIO_def.hpp"

