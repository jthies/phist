/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "TestWithType.h"

#ifdef PHIST_HAVE_SP

# include "phist_gen_s.h"
bool TestWithType<_ST_>::typeImplemented_ = false;

#ifdef PHIST_HAVE_CMPLX
# include "phist_gen_c.h"
bool TestWithType<_ST_>::typeImplemented_ = false;
#endif
#endif

# include "phist_gen_d.h"
bool TestWithType<_ST_>::typeImplemented_ = false;

#ifdef PHIST_HAVE_CMPLX
# include "phist_gen_z.h"
bool TestWithType<_ST_>::typeImplemented_ = false;
#endif
