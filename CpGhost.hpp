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
		printf("Before update: \n");
		ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);
		printf("After update: \n");
		return;
	}

	void write( const std::string * filename){
		printf("i will write now\n");
		asynData->toFile(asynData, (char *) (*filename).c_str(), 0); 
		return;
	}
	
	void read(const std::string * filename){
		printf("i will read now\n");
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
		//printf("Before update: \n");
	//	ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);
		//printf("After update: \n");
		return;
	}

	void write( const std::string * filename){
		printf("CpGhostDenseMatArray write function is not implemented: \n");
		//printf("i will write now\n");
	//	asynData->toFile(asynData, (char *) (*filename).c_str(), 0); 
		return;
	}
	
	void read(const std::string * filename){
		printf("CpGhostDenseMatArray read function is not implemented: \n");
		//printf("i will read now\n");
	//	asynData->fromFile(asynData, (char *) (*filename).c_str(), 0); 
	//	ghost_densemat_init_densemat(dataPtr, asynData, 0, 0);
		return;
	}
};




#endif
