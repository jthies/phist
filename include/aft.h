#ifndef __AFT_H__
#define __AFT_H__

#include "mpi.h"
#include "mpi-ext.h"
#include <unistd.h>
//#include "config.h"

// int verbose= 1;
int app_needs_repair(MPI_Comm *comm, char ** argv);
void errhandler_respawn(MPI_Comm* pcomm, int* errcode, ...);
int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm, char** argv);
void makeHostList(const char * machinefile, const int failCount, const int nd, const int nc, char * spawnHosts);

//int MPIX_Comm_replace(MPI_Comm worldwspares, MPI_Comm comm, MPI_Comm *newcomm);

#endif

// TODO: kill_all_procs_on_failed_processhosts
// active_machine_list -> kill_all_procs_on_failed_processhost
//
