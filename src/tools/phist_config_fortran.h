/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#ifndef PHIST_CONFIG_FORTRAN_H
#define PHIST_CONFIG_FORTRAN_H

/*! \file phist_config_fortran.h \ingroup tools

    This file contains additional macros that can be used in Fortran code to e.g. define kind parameters.
*/
#include "phist_config.h"
#ifdef PHIST_HAVE_GHOST
#include <ghost/config.h>
#endif

#ifndef PHIST_HAVE_GHOST
/*! \def G_LIDX_T 

kind parameter for ghost local indices (used to define row function arguments)
*/
#define G_LIDX_T C_INT32_T
# ifdef PHIST_FORCE_32BIT_GIDX
/*! \def G_GIDX_T 

kind parameter for ghost global indices (used to define row function arguments)
*/
# define G_GIDX_T C_INT32_T
# else
# define G_GIDX_T C_INT64_T
# endif
#else
# ifdef GHOST_IDX64_LOCAL
# define G_LIDX_T C_INT64_T
# else
# define G_LIDX_T C_INT32_T
# endif
# ifdef GHOST_IDX64_GLOBAL
# define G_GIDX_T C_INT64_T
# else
# define G_GIDX_T C_INT32_T
# endif
#endif

#if defined(PHIST_TESTING) && PHIST_OUTLEV >= 6
#define F_DEBUG 1
#endif

#endif
