/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_tools.h"
#include "phist_kernels.h"
#include "phist_core.h"
#include "phist_ScalarTraits.hpp"

#include "phist_iter_op.h"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_iter_op_def.hpp"
#include "phist_gen_c.h"
#include "phist_iter_op_def.hpp"
#endif

#include "phist_gen_d.h"
#include "phist_iter_op_def.hpp"
#include "phist_gen_z.h"
#include "phist_iter_op_def.hpp"

