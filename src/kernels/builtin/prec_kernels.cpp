/*! \file builtin/prec_kernels.cpp
 * wraps implementation of builtin kernel routines
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
#error "builtin kernels only work with MPI"
#endif

#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "phist_macros.h"
#include "phist_kernel_perfmodels.hpp"
#ifdef PHIST_HAVE_TEUCHOS
#include "phist_trilinos_macros.h"
#endif
#include "../phist_kernels.h"

#include "phist_typedefs.h"
#include "typedefs.hpp"
#include "phist_ScalarTraits.hpp"
extern "C" {
#include "bench_kernels.h"
}

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <cstring>
#include <sys/resource.h>

extern "C" {

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "../kernels_noprec.c"

#include "phist_gen_c.h"
#include "../kernels_noprec.c"
#endif

#include "phist_gen_z.h"
#include "../kernels_noprec.c"

#ifdef PHIST_HAVE_HIGH_PRECISION_KERNELS
#include "prec_kernels_d.h"
#endif


#include "phist_gen_d.h"
#ifdef PHIST_HAVE_HIGH_PRECISION_KERNELS
} //extern "C"
#include "prec_kernels_def.hpp"
#else
#include "../kernels_noprec.c"
} //extern "C"
#endif
