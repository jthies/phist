#ifndef __CPHELPERFUNCS_HPP__
#define __CPHELPERFUNCS_HPP__

#include <string>
#include <iostream>
#include <sstream>
#include <stdarg.h> 

#include <mpi.h>

#include "craftConf.h"


template <typename T>
std::string numberToString ( T number );

template <typename T>
T stringToNumber ( const std::string &text );

inline int mpiCommWorldRank();

inline int mpiCommWorldNumProcs();

void craftErr(const char *fmt, ...);
void craftDbg(int level, const char *fmt, ...);
void craftAbort(int rc, const char *fmt, ...);
void craftTime(const std::string str);
void craftTime(const std::string str, const MPI_Comm * const comm);
void getEnvVal(int &var, const char * str);
void checkDirectoryName(std::string* const s);
void checkPathName(std::string* const s);
int craftLog(MPI_Comm const * comm, const char *fmt, ...);

#endif

