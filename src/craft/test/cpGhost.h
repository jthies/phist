#ifndef __CPGHOST_H__
#define __CPGHOST_H__

#include "cpDataBase.h"
#include "ghost.h"

class CpGhostDenseMat: public CpDataBase  
{
private:

	ghost_densemat * dataPtr;
	ghost_densemat * asynData;

public:
	CpGhostDenseMat();
	
	CpGhostDenseMat( ghost_densemat * asynCopy_);

	~CpGhostDenseMat();	

	int update(){
		printf("Before update: \n");
		ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);
		printf("After update: \n");
		return 0;
	}

	int write(char * filename){
		printf("i will write now\n");
		asynData->toFile(asynData, filename, 0); 
		return 0;
	}
	
	int read(char * filename){
		printf("i will read now\n");
		asynData->fromFile(asynData, filename, 0); 
		ghost_densemat_init_densemat(dataPtr, asynData, 0, 0);
		return 0;
	}
};


CpGhostDenseMat::CpGhostDenseMat(){

}

CpGhostDenseMat::CpGhostDenseMat( ghost_densemat * asynCopy_){
	asynData = new ghost_densemat[1];
	dataPtr  = new ghost_densemat[1];
	
	dataPtr = asynCopy_;
	ghost_densemat_create(&asynData, dataPtr->context, dataPtr->traits); 
	ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);

//	asynData->copy(dataPtr, asynData);											// (*source, *destination)	
}


CpGhostDenseMat::~CpGhostDenseMat(){

}


#endif
