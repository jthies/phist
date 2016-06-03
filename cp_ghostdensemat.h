/**
 * Author: Faisal Shahzad
 * This file contains two classes for denseMat and denseMatArrays.
 *
 */

#ifndef __CP_GHOSTDENSEMAT_H__
#define __CP_GHOSTDENSEMAT_H__

#include "enum.h"

class CpGhostDenseMat 
{
public:
		CpGhostDenseMat();
		~CpGhostDenseMat();
		int add(const std::string name, ghost_densemat * const item, const MPI_Comm FtComm, const std::string cpPath_ );
		int update();	// check it the item is defined. if it is defined, update according to nRows and array
	   	int write();	
		int read();
private:
		std::string cpPath;
		MPI_Comm cpMpiComm;
		std::string name;
		ghost_densemat * denseMatAsync;			// the poointer of acsynchronous vector to be written
	 	ghost_densemat * denseMat;		// the pointer of actual vetor 
};

CpGhostDenseMat::CpGhostDenseMat()
{
	cpMpiComm = MPI_COMM_WORLD;
	cpPath = "";
	name = "";
}

CpGhostDenseMat::~CpGhostDenseMat()
{

}

int CpGhostDenseMat::add(const std::string name_, ghost_densemat * const item, const MPI_Comm FtComm, const std::string cpPath_)
{
	denseMat = item;	
	name = name_;	
	cpMpiComm = FtComm;
	cpPath = cpPath_;
	denseMatAsync = new ghost_densemat[1];
	ghost_densemat_create(&denseMatAsync, item->context, item->traits);
	ghost_densemat_init_densemat(denseMatAsync, item, 0, 0);

	return 0;
}

int CpGhostDenseMat::update(){
	std::cout << "updating DenseMat:" << name << std::endl;
	ghost_densemat_init_densemat(denseMatAsync, denseMat, 0, 0);
//			char * vecstr;
//			ghost_densemat_string(&vecstr, itAsync->second);
//			printf("y: \n%s\n", vecstr);	
	return 0;
}

int CpGhostDenseMat::write(){
	std::cout << "writing file now " << name << std::endl ;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);		
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);
	denseMatAsync->toFile(denseMatAsync, filename, 0);
	std::cout << "writing done " << name << std::endl ;
	return 0;
}

int CpGhostDenseMat::read(){
	std::cout << "read file now" << name << std::endl;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);
	denseMatAsync->fromFile(denseMatAsync, filename, 0);
	ghost_densemat_init_densemat(denseMat, denseMatAsync, 0, 0);
	char * readstr = new char[256];	
	ghost_densemat_string (&readstr, denseMat);
	
	return 0;
}

// ===== GHOST DENSE-MAT-ARRAY CLASS ===== //
//
class CpGhostDenseMatArray
{
public:
		CpGhostDenseMatArray();
		~CpGhostDenseMatArray();
		int add(const std::string name, ghost_densemat ** const item, const size_t nDenseMat_, const MPI_Comm FtComm, const std::string cpPath_, const int toCpDenseMat_);
		int update();	// check it the item is defined. if it is defined, update according to nRows and array
	   	int write();	
		int read();
private:
		std::string cpPath;
		MPI_Comm cpMpiComm;
		std::string name;
		ghost_densemat ** denseMatArrayAsync;			// the poointer of acsynchronous vector to be written
	 	ghost_densemat ** denseMatArray;		// the pointer of actual vetor 
		size_t nDenseMat;
		int toCpDenseMat;
};

CpGhostDenseMatArray::CpGhostDenseMatArray()
{
	cpMpiComm = MPI_COMM_WORLD;
	cpPath = "";
	name = "";
	nDenseMat = 0;
	toCpDenseMat = 0;
}

CpGhostDenseMatArray::~CpGhostDenseMatArray()
{

}

int CpGhostDenseMatArray::add(const std::string name_, ghost_densemat ** const item, const size_t nDenseMat_, const MPI_Comm FtComm, const std::string cpPath_, const int toCpDenseMat_)
{
	denseMatArray = item;
	nDenseMat = nDenseMat_;	
	name = name_;	
	cpMpiComm = FtComm;
	cpPath = cpPath_;
	denseMatArrayAsync = new ghost_densemat*[nDenseMat];
	toCpDenseMat = toCpDenseMat_;
	for(size_t i = 0; i < nDenseMat ; ++i)
	{
		ghost_densemat_create( &denseMatArrayAsync[i], item[i]->context, item[i]->traits);
		ghost_densemat_init_densemat ( denseMatArrayAsync[i], denseMatArray[i], 0, 0);	
	}
	
	return 0;
}

int CpGhostDenseMatArray::update(){
	std::cout << "updating DenseMatArray: " << name << std::endl;
	
	for(size_t i = 0; i < nDenseMat ; ++i)
	{
		ghost_densemat_init_densemat ( denseMatArrayAsync[i], denseMatArray[i], 0, 0);	
	}	
	/*
	ghost_densemat_init_densemat(denseMatAsync, denseMat, 0, 0);
	std::cout << "Densemat is updated"  << std::endl;
	*/
//			char * vecstr;
//			ghost_densemat_string(&vecstr, itAsync->second);
//			printf("y: \n%s\n", vecstr);	
	return 0;
}

int CpGhostDenseMatArray::write(){
	std::cout << "writing DenseMatArray: " << name << std::endl ;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);	
	
	if(toCpDenseMat == ALL){
		std::cout << "writing all MULTIVEC" << std::endl;
		for(int i = 0; i < nDenseMat ; ++i){
			char * filename = new char[256];
			sprintf(filename, "%s/%s%d-rank%d.cp", cpPath.c_str(), name.c_str(), i, myrank);

			FILE * fp;
			if( NULL == (fp = fopen(filename, "w+")) ) {
				fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
				return -1;
			}		
			denseMatArrayAsync[i]->toFile(denseMatArrayAsync[i], filename, 0);
		}
	}
	if(toCpDenseMat == CYCLIC){
		static size_t cyclicCpCounter = 0;	//  which Densemat is to be checkpoinited in the cyclic case. 	
		std::cout << "writing all MULTIVEC" << std::endl;
		char * filename = new char[256];
		sprintf(filename, "%s/%s%zu-rank%d.cp", cpPath.c_str(), name.c_str(), cyclicCpCounter, myrank);

		FILE * fp;
		if( NULL == (fp = fopen(filename, "w+")) ) {
			fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
			return -1;
		}		
		denseMatArrayAsync[cyclicCpCounter]->toFile(denseMatArrayAsync[cyclicCpCounter], filename, 0);
		cyclicCpCounter++;	
		if(cyclicCpCounter == nDenseMat)
		{
			cyclicCpCounter = 0;
		}
	}
	if( toCpDenseMat >= 0){
		std::cout << "writing toCpDenseMat" << toCpDenseMat<< std::endl;
		char * filename = new char[256];
		sprintf(filename, "%s/%s%d-rank%d.cp", cpPath.c_str(), name.c_str(), toCpDenseMat, myrank);

		FILE * fp;
		if( NULL == (fp = fopen(filename, "w+")) ) {
			fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
			return -1;
		}		
		denseMatArrayAsync[toCpDenseMat]->toFile(denseMatArrayAsync[toCpDenseMat], filename, 0);
	}
	return 0;
}

int CpGhostDenseMatArray::read(){
	std::cout << "read file now" << name << std::endl;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);	
	
	if(toCpDenseMat == ALL){
		std::cout << "reading all MULTIVEC" << std::endl;
		for(int i = 0; i < nDenseMat ; ++i){
			char * filename = new char[256];
			sprintf(filename, "%s/%s%d-rank%d.cp", cpPath.c_str(), name.c_str(), i, myrank);

			FILE * fp;
			if( NULL == (fp = fopen(filename, "r")) ) {
				fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
				return -1;
			}		
			denseMatArrayAsync[i]->fromFile(denseMatArrayAsync[i], filename, 0);
			ghost_densemat_init_densemat ( denseMatArray[i], denseMatArrayAsync[i], 0, 0);	
		}
	}
	if(toCpDenseMat == CYCLIC){
		// TODO: 	check the restart with cyclic case. restart will have to determine which
		std::cerr << "Restart with CYCLIC case is not reading yet." << std::endl;
		/*static size_t cyclicCpCounter = 0;	//  which Densemat is to be checkpoinited in the cyclic case. 	
		std::cout << "writing all MULTIVEC" << std::endl;
		char * filename = new char[256];
		sprintf(filename, "%s/%s%zu-rank%d.cp", cpPath.c_str(), name.c_str(), cyclicCpCounter, myrank);

		FILE * fp;
		if( NULL == (fp = fopen(filename, "r")) ) {
			fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
			return -1;
		}		
		denseMatArrayAsync[cyclicCpCounter]->fromFile(denseMatArrayAsync[cyclicCpCounter], filename, 0);
		cyclicCpCounter++;	
		if(cyclicCpCounter == nDenseMat)
		{
			cyclicCpCounter = 0;
		}*/
	}
	if( toCpDenseMat >= 0){
		std::cout << "writing toCpDenseMat" << toCpDenseMat<< std::endl;
		char * filename = new char[256];
		sprintf(filename, "%s/%s%d-rank%d.cp", cpPath.c_str(), name.c_str(), toCpDenseMat, myrank);

		FILE * fp;
		if( NULL == (fp = fopen(filename, "r")) ) {
			fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
			return -1;
		}		
		denseMatArrayAsync[toCpDenseMat]->fromFile(denseMatArrayAsync[toCpDenseMat], filename, 0);
		ghost_densemat_init_densemat ( denseMatArray[toCpDenseMat], denseMatArrayAsync[toCpDenseMat], 0, 0);	
	}
	return 0;
}



#endif
