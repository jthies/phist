#ifndef __CHECKPOINT_HPP__
#define __CHECKPOINT_HPP__

#include <string>

#include "enum.h"
#include "CpBase.hpp"
#include "CpPOD.hpp"
#include "CpGhost.hpp"
#include "CpArray.hpp"

#include <string>
#include <sstream>
#include <map>
#include <mpi.h>

#ifdef SCR
extern "C"{
	#include <scr.h>
}
#endif

class Checkpoint
{
protected:
	MPI_Comm cpMpiComm;
	std::string cpPath;
	bool useSCR;
	bool cpCommitted;

public:
	Checkpoint();  
	~Checkpoint();  
  	void setCpPath (const std::string cpPath_);
	void setComm (const MPI_Comm cpMpiComm_);
	void enableSCR();
	void commit();

	typedef std::map<const std::string,CpBase *> cp_const_map;		// TODO: check if they can be privatly declared
  	typedef std::map<const std::string,CpBase *> cp_copy_map;
 	// contains points to user workspace
  	cp_const_map objects;
  	// contains asynchronous copies of objects
  	cp_copy_map async_copies;

 
  	void read();	
  	void write();
	void update();
  
// implementation of add() for anything that is CpBase,
// anything else will give an error message
void add(std::string label, CpBase * p)
{
    objects[label] = p;
}

// specialized implementation of add for POD (plain old data), which is obviously
// not derived from our base class CpBase. The same can be done for e.g. 
// ghost_densemat
// TODO: here a new object is created that will not be deleted unless we have some
//       memory management in the background, the "normal" add function just adds 
//       the pointer to the map and the person who created the object is responsible
//       for deleting it.

	// ===== POD ===== //
  	void add(std::string label, int * const i){this->add(label, new CpPOD<int>(i));}
  	void add(std::string label, float * const i){this->add(label, new CpPOD<float>(i));}
	void add(std::string label, double * const d){this->add(label, new CpPOD<double>(d));}

	// ===== GHOST DENSE MATRIX ===== //
	void add(std::string label, ghost_densemat * const GDM)
	{
		this->add(label, new CpGhostDenseMat(GDM));
	}
	void add(std::string label, ghost_densemat ** const GDMArray, const size_t nDenseMat_, const int toCpDenseMat_=ALL)
	{
		this->add(label, new CpGhostDenseMatArray(GDMArray, nDenseMat_, toCpDenseMat_) );
	}

	// ===== POD ARRAY ===== // 
	template <class T>
	void add(std::string label, T* const arrayPtr_, const size_t nRows_){
			this->add(label, new CpArray<T>(arrayPtr_, nRows_));
	}

#endif

};


Checkpoint::Checkpoint(){
	cpMpiComm 	= MPI_COMM_WORLD;	
	cpPath 		= "";
	cpCommitted = false;
	useSCR 		= false;
}	

Checkpoint::~Checkpoint(){
#ifdef SCR
	if(useSCR == true){
		printf("====== SCR_Finalize is done ====== \n");
		SCR_Finalize();
	}
#endif
}

void Checkpoint::setCpPath(const std::string cpPath_){
	cpPath = cpPath_;
	return;
}
	
void Checkpoint::setComm(const MPI_Comm cpMpiComm_){
	cpMpiComm = cpMpiComm_;
#ifdef SCR
	if(useSCR == true){
		SCR_Init(&cpMpiComm);
	}
#endif
	return;
}

void Checkpoint::enableSCR(){
#ifdef SCR
	useSCR = true;
//	printf("useSCR is set true \n");

#else 
	printf("ERROR: Checkpoint-lib is not compiled with SCR. Checkpoints will be written at default path\n");
#endif
}

void Checkpoint::commit(){
	cpCommitted = true;	
}


