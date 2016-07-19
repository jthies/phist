#ifndef __CPPHISTMVEC_H__
#define __CPPHISTMVEC_H__

#include "CpBase.hpp"
#include "ghost.h"

#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_enums.h"
#include "phist_macros.h"
#include "phist_lapack.h"
//#include "phist_operator_decl.h"

//#include "phist_gen_z.h"

 #include "phist_ScalarTraits.hpp"
typedef phist::ScalarTraits<_ST_> st;


class CpPhistMvec: public CpBase  
{
private:

	TYPE(mvec_ptr) dataPtr;
	TYPE(mvec_ptr) asynData;
	int iflag;

public:
	
	CpPhistMvec(TYPE(mvec_ptr)  dataPtr_){
		printf("CpPhistMvec constructor is called here: \n");
		dataPtr = dataPtr_;
		int nVec;
		iflag = 0;
		PHIST_CHK_IERR(SUBR(mvec_num_vectors)(dataPtr,&nVec,&iflag),iflag);
		printf("vectors in dataPtr: %d\n",nVec);
		phist_const_map_ptr map;
		PHIST_CHK_IERR(SUBR(mvec_get_map)(dataPtr, &map, &iflag), iflag); 
  		PHIST_CHK_IERR(SUBR(mvec_create)(&asynData,map,nVec,&iflag),iflag);


  	//	PHIST_CHK_IERR(SUBR(mvec_create)(&asynData,A_op->domain_map,1,iflag),*iflag);

		//asynData = new ghost_densemat[1];
		//dataPtr  = new ghost_densemat[1];
	
		//ghost_densemat_create(&asynData, dataPtr->context, dataPtr->traits); 
		//ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);


	}
	~CpPhistMvec(){}	

	void update(){
		printf("CpPhistMvec Before update: \n");
		PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),dataPtr,st::zero(),asynData,&iflag),iflag);
		//ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);
		//printf("After update: \n");
		return;
	}

	void write( const std::string * filename){
		printf("CpGhostDenseMatArray i will write now\n");
		PHIST_CHK_IERR(SUBR(mvec_write_bin)(asynData, *filename ,&iflag),iflag);
		//asynData->toFile(asynData, (char *) (*filename).c_str(), 0); 
		return;
	}
	
	void read(const std::string * filename){
		printf("i will read now\n");
		//asynData->fromFile(asynData, (char *) (*filename).c_str(), 0); 
		//ghost_densemat_init_densemat(dataPtr, asynData, 0, 0);
		return;
	}
};

/*
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

*/


#endif
