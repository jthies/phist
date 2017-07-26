/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_KPM_H
#define PHIST_KPM_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_void_aliases.h"
#include "phist_core_flags.h"

#endif //DOXYGEN

#ifdef __cplusplus
extern "C" {
#endif
#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_kpm_decl.h"
#include "phist_gen_c.h"
#include "phist_kpm_decl.h"
#endif
#include "phist_gen_d.h"
#include "phist_kpm_decl.h"
#include "phist_gen_z.h"
#include "phist_kpm_decl.h"

#include "phist_gen_clean.h"

#ifdef __cplusplus
}
#endif

#endif
