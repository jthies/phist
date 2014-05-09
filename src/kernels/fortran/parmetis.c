/*! \file parmetis.c
 * \brief defines the parmetis wrapper function parmetis_v3_partkway_f_
 * \author Achim Basermann, Melven Zoellner
 *
 * $Id: parmetis.c 89 2011-09-27 13:01:01Z zoel_me $
 *
 */

#include "phist_config.h"

#ifdef PHIST_HAVE_PARMETIS

#include <parmetis.h>

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
                            float * j,
                            float * k,
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
                            float * j,
                            float * k,
                            int64_t * l,
                            int64_t * m,
                            int64_t * n,
                            int * fComm)
{
  MPI_Comm cComm = MPI_Comm_f2c((MPI_Fint)*fComm);
  MPI_Comm *comm = &cComm;

  /* the RefineKway doesn't support the uncoupling the vertices from the processors int the current version, so we try AdaptiveRepart */
  /* newer versions fail in the above function for empty processes */
  float itr = 1000;
  ParMETIS_V3_AdaptiveRepart(a,b,c,d,NULL,e,f,g,h,i,j,k,&itr,l,m,n,comm);
  /*ParMETIS_V3_RefineKway(a,b,c,d,e,f,g,h,i,j,k,l,m,n,comm);*/
}

#endif

