#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_orthog.h"
#include "phist_kernels.h"

#include "phist_ScalarTraits.hpp"

//#ifdef PHIST_KERNEL_LIB_TPETRA
//#define USE_TRILINOS_ORTHO_MANAGER
//#endif

#ifdef USE_TRILINOS_ORTHO_MANAGER

#include "phist_rcp_helpers.hpp"
#include "phist_operator.h"
#include "phist_BelosOperatorTraits.hpp"

#ifdef PHIST_KERNEL_LIB_TPETRA
#include "BelosTpetraAdapter.hpp"
#include "phist_tpetra_typedefs.hpp"
#else
#error "not implemented"
#endif

#include "BelosOrthoManager.hpp"
// we can try any of the methods implemented in Belos:
#include "BelosTsqrOrthoManager.hpp"
#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "trili_orthog_def.hpp"
#include "phist_gen_c.h"
#include "trili_orthog_def.hpp"
#endif
#include "phist_gen_d.h"
#include "trili_orthog_def.hpp"
#include "phist_gen_z.h"
#include "trili_orthog_def.hpp"
#else
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_orthog_def.hpp"
#include "phist_gen_c.h"
#include "phist_orthog_def.hpp"
#endif
#include "phist_gen_d.h"
#include "phist_orthog_def.hpp"
#include "phist_gen_z.h"
#include "phist_orthog_def.hpp"
#endif
