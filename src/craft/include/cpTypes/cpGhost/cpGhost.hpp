#ifndef __CPGHOST_HPP__
#define __CPGHOST_HPP__

#include "../../cpEnum.h"
#include "../../cpBase.hpp"
#include "ghost.h"


class CpGhostDenseMat: public CpBase  
{
public:
	CpGhostDenseMat(ghost_densemat *  dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpGhostDenseMat(){}	
	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:
	ghost_densemat *  dataPtr;
	ghost_densemat * asynData;
};



class CpGhostDenseMatArray: public CpBase  
{
public:
	CpGhostDenseMatArray(ghost_densemat **  dataPtr_, const size_t nDenseMat_, const int toCpDenseMat_=ALL, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpGhostDenseMatArray(){}	

	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:
	ghost_densemat ** dataPtr;
	ghost_densemat ** asynData;
	size_t nDenseMat;
	int toCpDenseMat; 
	size_t cyclicCpCounter;
  MPI_Comm cpMpiComm ;
};



class CpGhostSparseMat: public CpBase  	// TODO: not fully implemented
{
public:
	CpGhostSparseMat(ghost_sparsemat *  dataPtr_);
	~CpGhostSparseMat(){}	
	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:;
	ghost_sparsemat * dataPtr;
	ghost_sparsemat * asynData;
};


#endif
