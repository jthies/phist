#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"
#include "phist_precon.h"
#include "phist_macros.h"
#include "phist_ScalarTraits.hpp"
#include <stdlib.h>

#include "phist_PreconTraits.hpp"

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_precon_def.hpp"
#include "phist_gen_c.h"
#include "phist_precon_def.hpp"
#endif
#include "phist_gen_d.h"
#ifdef PHIST_KERNEL_LIB_EPETRA
# ifdef PHIST_HAVE_IFPACK
# include "tpl/phist_Ifpack_def.hpp"
# endif
#endif
#include "phist_precon_def.hpp"
#include "phist_gen_z.h"
#include "phist_precon_def.hpp"
