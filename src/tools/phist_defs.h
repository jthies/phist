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

/*! \file phist_defs.h \ingroup tools

this file contains general macros for common return values
and output levels etc. and can be included in C, C++ or Fortran 
code alike */

/*! \name flags for output \ingroup tools
 *  these can be passed to PHIST_OUT(FLAG,...)
 *  to get a certain amount of coherence in the 
 *  way screen output is handled. The PHIST_OUTLEV
 *  macro defines which messages are printed and
 */
/*!@{*/
 
/*! \def PHIST_ERROR */
#define PHIST_ERROR   0
/*! \def PHIST_WARNING */
#define PHIST_WARNING 1
/*! \def PHIST_INFO */
#define PHIST_INFO 2
/*! \def PHIST_VERBOSE */
#define PHIST_VERBOSE 3
/*! \def PHIST_EXTREME */
#define PHIST_EXTREME 4
/*! \def PHIST_DEBUG */
#define PHIST_DEBUG 5
/*! \def PHIST_TRACE */
#define PHIST_TRACE 6
/*!@}*/
#if defined (PHIST_PERFCHECK)
#define PHIST_PERFWARNING PHIST_WARNING
#else
#define PHIST_PERFWARNING PHIST_VERBOSE
#endif

/*! \name return types \ingroup tools
 * error messages
 */
/*!@{*/
    
/*! \def PHIST_SUCCESS */
#define PHIST_SUCCESS 0
/*! \def PHIST_FUNCTIONAL_ERROR */
#define PHIST_FUNCTIONAL_ERROR -1
/*! \def PHIST_MPI_ERROR */
#define PHIST_MPI_ERROR -33
/*! \def PHIST_MEM_ALLOC_FAILED */
#define PHIST_MEM_ALLOC_FAILED -44
/*! \def PHIST_INVALID_INPUT */
#define PHIST_INVALID_INPUT -55
/*! \def PHIST_INTEGER_OVERFLOW */
#define PHIST_INTEGER_OVERFLOW -66
/*! \def PHIST_CAUGHT_EXCEPTION */
#define PHIST_CAUGHT_EXCEPTION -77
/*! \def PHIST_BAD_CAST */
#define PHIST_BAD_CAST -88
/*! \def PHIST_NOT_IMPLEMENTED */
#define PHIST_NOT_IMPLEMENTED -99
/*! \def PHIST_DEPRECATED */
#define PHIST_DEPRECATED +99
/*!@}*/

#endif
