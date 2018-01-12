/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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

# ifdef PHIST_HAVE_CMPLX
# include "phist_gen_c.h"
# include "MatrixIO_decl.h"
# endif

#endif

#include "phist_gen_d.h"
#include "MatrixIO_decl.h"

#ifdef PHIST_HAVE_CMPLX

#include "phist_gen_z.h"
#include "MatrixIO_decl.h"

#endif
  }
}

#endif
