#ifndef __CPARRAY_HPP__
#define __CPARRAY_HPP__

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "cpHelperFuncs.hpp"
#include "cpEnum.h"
#include "cpBase.hpp"



// ===== POD ARRAY ===== // 
template <class T>
class CpArray: public CpBase
{
public:
  CpArray( T * const dataPtr_, const size_t nRows_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpArray(){}	

	int read(const std::string * filename);
  int write(const std::string * filename);
	int update();
	void print();
	
private:
	T * dataPtr;
	T * asynData;
	size_t nRows;
	void copyArray(T * const src, T * const dst, const size_t numElems);

  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
	
};



// ===== POD MULTI-ARRAY ===== // 
template <class T>
class CpMultiArray: public CpBase
{
public:
	CpMultiArray( T ** const dataPtr_, const size_t nRows_, const size_t nCols_, 
                const int toCpCol_=ALL, 
                const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpMultiArray(){}	

	int read(const std::string * filename);
	int write(const std::string * filename);
  int update();
	void print();

private:
	T ** dataPtr;
	T ** asynData;
	size_t nRows;
	size_t nCols;
	int toCpCol;
	size_t cyclicCpCounter;
	void copyMultiArray(T ** const src, T ** const dst, const size_t nRows, const size_t nCols);

	int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);

};

#endif

