/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file parmetis.c
 * Provides a wrapper for Fortran for the C function parmetis_v3_partkway from ParMetis
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_HAVE_PARMETIS

#include <parmetis.h>
#include <stdlib.h>

#if REALTYPEWIDTH == 32
typedef float real;
#elif REALTYPEWITH == 64
typedef doublle real;
#else
#error "do not find or understand RELTYPEWIDTH, should be defined in metis.h as 32 or 64"
#endif

/*! wrapper for parmetis, so it can be called from fortran:<br>
 * it is necessary to convert the mpi communicator from c
 * to fortran (and this is only possible in c!)<br>
 * description of parameters: see ParMETIS reference
 */
void parmetis_v3_partkway_f_ (int64_t * a,
                            int64_t * b,
                            int64_t * c,
                            int64_t * d,
                            int64_t * e,
                            int64_t * f,
                            int64_t * g,
                            int64_t * h,
                            int64_t * i,
                            real * j,
                            real * k,
                            int64_t * l,
                            int64_t * m,
                            int64_t * n,
                            int * fComm)
{
  MPI_Comm cComm = MPI_Comm_f2c((MPI_Fint)*fComm);
  MPI_Comm *comm = &cComm;

 ParMETIS_V3_PartKway(a,b,c,d,e,f,g,h,i,j,k,l,m,n,comm);
 /*ParMETIS_PartKway(a,b,c,d,e,f,g,i,l,m,n,comm);*/
}

void parmetis_v3_refinekway_f_ (int64_t * a,
                            int64_t * b,
                            int64_t * c,
                            int64_t * d,
                            int64_t * e,
                            int64_t * f,
                            int64_t * g,
                            int64_t * h,
                            int64_t * i,
                            real * j,
                            real * k,
                            int64_t * l,
                            int64_t * m,
                            int64_t * n,
                            int * fComm)
{
  MPI_Comm cComm = MPI_Comm_f2c((MPI_Fint)*fComm);
  MPI_Comm *comm = &cComm;

  /* the RefineKway doesn't support the uncoupling the vertices from the processors int the current version, so we try AdaptiveRepart */
  /* newer versions fail in the above function for empty processes */
  real itr = 1000;
  ParMETIS_V3_AdaptiveRepart(a,b,c,d,NULL,e,f,g,h,i,j,k,&itr,l,m,n,comm);
  /*ParMETIS_V3_RefineKway(a,b,c,d,e,f,g,h,i,j,k,l,m,n,comm);*/
}

#endif

