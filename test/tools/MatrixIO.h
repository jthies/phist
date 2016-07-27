#ifndef PHIST_TEST_MATRIXIO_H
#define PHIST_TEST_MATRIXIO_H

#include "phist_kernels.h"

namespace phist
{
  namespace testing
  {

    // set the flag for reading/creating the matrix dpeending on #rows and #vectors
    int getSparseMatCreateFlag(int N, int NV);

#ifdef PHIST_HAVE_SP

#include "phist_gen_s.h"
#include "MatrixIO_decl.h"

#include "phist_gen_c.h"
#include "MatrixIO_decl.h"

#endif

#include "phist_gen_d.h"
#include "MatrixIO_decl.h"

#include "phist_gen_z.h"
#include "MatrixIO_decl.h"

  }
}

#endif
