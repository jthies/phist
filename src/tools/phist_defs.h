/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_DEFS_H
#define PHIST_DEFS_H

/* this file contains general macros for common return values
and output levels etc. and can be included in C, C++ or Fortran 
code alike */

/* these can be passed to PHIST_OUT(FLAG,...)
 to get a certain amount of coherence in the 
 way screen output is handled. The PHIST_OUTLEV
 macro defines which messages are printed and
 which aren't. */
 
#define PHIST_ERROR   0
#define PHIST_WARNING 1
#define PHIST_INFO 2
#define PHIST_VERBOSE 3
#define PHIST_EXTREME 4
#define PHIST_DEBUG 5
#define PHIST_TRACE 6

#if defined (PHIST_PERFCHECK)
#define PHIST_PERFWARNING PHIST_WARNING
#else
#define PHIST_PERFWARNING PHIST_INFO
#endif

/* return types */
#define PHIST_SUCCESS 0
#define PHIST_FUNCTIONAL_ERROR -1
#define PHIST_MPI_ERROR -33
#define PHIST_MEM_ALLOC_FAILED -44
#define PHIST_INVALID_INPUT -55
#define PHIST_INTEGER_OVERFLOW -66
#define PHIST_CAUGHT_EXCEPTION -77
#define PHIST_BAD_CAST -88
#define PHIST_NOT_IMPLEMENTED -99
#define PHIST_DEPRECATED +99

#endif
