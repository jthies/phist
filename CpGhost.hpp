#ifndef __CPGHOST_H__
#define __CPGHOST_H__

#include "CpBase.hpp"
#include "ghost.h"

class CpGhostDenseMat: public CpBase  
{
private:

	ghost_densemat *  dataPtr;
	ghost_densemat * asynData;

public:
	CpGhostDenseMat(ghost_densemat *  dataPtr_){
		asynData = new ghost_densemat[1];
		dataPtr  = new ghost_densemat[1];
	
		dataPtr = dataPtr_;
		ghost_densemat_create(&asynData, dataPtr->context, dataPtr->traits); 
		ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);


	}
	~CpGhostDenseMat(){}	

	void update(){
		printf("Updating ghost densemat: \n");
		ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);
		return;
	}

	void write( const std::string * filename){
		printf("Writing ghost densemat: \n");
		asynData->toFile(asynData, (char *) (*filename).c_str(), MPI_COMM_SELF);	// TODO:  enbale local and private writes.
		return;
	}
	
	void read(const std::string * filename){
		printf("Reading ghost densemat: \n");
		asynData->fromFile(asynData, (char *) (*filename).c_str(), 0); 
		ghost_densemat_init_densemat(dataPtr, asynData, 0, 0);
		return;
	}
};

class CpGhostDenseMatArray: public CpBase  
{
private:

	ghost_densemat ** dataPtr;
	ghost_densemat ** asynData;
	size_t nDenseMat;
	int toCpDenseMat; 
public:
	CpGhostDenseMatArray(ghost_densemat **  dataPtr_, const size_t nDenseMat_, const int toCpDenseMat_){
		dataPtr = dataPtr_;
		nDenseMat = nDenseMat_;
		toCpDenseMat = toCpDenseMat_;
		
		asynData = new ghost_densemat*[nDenseMat];
		for(size_t i = 0; i < nDenseMat ; ++i)
		{
			ghost_densemat_create( &asynData[i], dataPtr[i]->context, dataPtr[i]->traits);
			ghost_densemat_init_densemat ( asynData[i], dataPtr[i], 0, 0);	
		}
	}

	~CpGhostDenseMatArray(){}	

	void update(){
		printf("CpGhostDenseMatArray update function is not implemented: \n");
		for(size_t i = 0; i < nDenseMat ; ++i)
		{
			printf("Updating ghost densemat: %d\n", i);
			ghost_densemat_init_densemat ( asynData[i], dataPtr[i], 0, 0);	
		}
		return;
	}

	void write( const std::string * filename){
		//printf("i will write now\n");
		for(size_t i = 0; i < nDenseMat ; ++i)
		{
			printf("Writing ghost densemat: %d\n", i);
			asynData[i]->toFile(asynData[i], (char *) (*filename).c_str(), MPI_COMM_SELF); // TODO:  enbale local and private writes.
		}
		return;
	}
	
	void read(const std::string * filename){
		printf("CpGhostDenseMatArray read function is not implemented: \n");
		for(size_t i = 0; i < nDenseMat ; ++i)
		{
			printf("Writing ghost densemat: %d\n", i);
			asynData[i]->fromFile(asynData[i], (char *) (*filename).c_str(), 0); 
			ghost_densemat_init_densemat(dataPtr[i], asynData[i], 0, 0);
		}
		return;
	}
};




#endif
