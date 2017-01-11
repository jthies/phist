#ifndef __AFT_H__
#define __AFT_H__

#include <unistd.h>
#include <string>
#include <vector>
#include <set>

#include "mpi.h"
#include "mpi-ext.h"

#include "craftConf.h"
#include "cpHelperFuncs.hpp"
#include "cpEnum.h"

int app_needs_repair(MPI_Comm *comm, char ** argv);
void errhandler_respawn(MPI_Comm* pcomm, int* errcode, ...);
int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm, char** argv);   // NONSHRINKING recovery ( REUSE, NOREUSE node policy)

int initRescueNodeList(MPI_Comm * const comm);
int writeActiveMachineList(std::vector<std::string> *activeNodeList_, MPI_Comm * const comm);
int makeActiveMachineList(std::vector<std::string> *activeNodeList_, MPI_Comm * const comm);

int makeSpawnList(const int nd, const int * const failedRanks, MPI_Comm * const comm, std::string * spawnList);
int makeSpawnListReuse(const int nd, const int * const failedRanks, MPI_Comm * const comm, std::string * spawnList);
int makeSpawnListNoReuse(const int nd, const int * const failedRanks, MPI_Comm * const comm, std::string * spawnList);

int writeFailedList(const int nd, const int * const failedRanks, MPI_Comm * const comm);

int getNumFailed(const MPI_Comm * const comm, const MPI_Comm * scomm);
int getFailedRanks(const MPI_Comm * const comm, const MPI_Comm * scomm, int * failedRanks);

int removeMachineFiles(MPI_Comm * const comm);
int removeFile(const char * filename, MPI_Comm * const comm);

int printNodeName( const MPI_Comm * const comm);
int getFirstRank(MPI_Comm* const comm);
int copyVecToSet(std::vector<std::string> * vecSrc, std::set<std::string> *setDst);


int findSetDiffence(std::set<std::string> * s1, std::set<std::string> * s2, std::set<std::string> * res);
int printSet(std::string toPrint, std::set<std::string> s);
int writeSetList(const std::string setFileName, std::set<std::string> list, MPI_Comm * const comm);
int writeVectorList(const std::string vecFileName, std::vector<std::string> vec, MPI_Comm * const comm);

const char * machinefileActiveProcs = "machinefileActiveProcs";
const char * machinefileFailedProcs = "machinefileFailedProcs";
const char * machinefileSpawnProcs  = "machinefileSpawnProcs";
const char * machinefileRescueProcs  = "machinefileRescueProcs";


#endif

// TODO: kill_all_procs_on_failed_processhosts
// active_machine_list -> kill_all_procs_on_failed_processhost
//
