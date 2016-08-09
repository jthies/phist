#ifndef PHIST_SIMPLE_LANCZOS_FT_H
#define PHIST_SIMPLE_LANCZOS_FT_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_operator.h"
#include "phist_enums.h"

#include "cp_options.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif	// CPP
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_simple_lanczos_decl_ft.h"
#include "phist_gen_c.h"
#include "phist_simple_lanczos_decl_ft.h"
#endif
#include "phist_gen_d.h"
#include "phist_simple_lanczos_decl_ft.h"
#include "phist_gen_z.h"
#include "phist_simple_lanczos_decl_ft.h"
#include "phist_gen_clean.h"
#ifdef __cplusplus
}
#endif	// CPP

#endif	// header
