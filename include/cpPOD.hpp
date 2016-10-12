#ifndef __CPPOD_HPP__
#define __CPPOD_HPP__

#include "cpHelperFuncs.hpp"
#include "cpBase.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

template <class T>
class CpPOD: public CpBase
{
public:
	CpPOD( T * const dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpPOD(){}	

	int read(const std::string * filename);
	int write(const std::string * filename);
	int update();
	
private:
	T * dataPtr;
	T * asynData;
  
  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
	
};




#endif
