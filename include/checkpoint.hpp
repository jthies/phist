#ifndef __CHECKPOINT_HPP__
#define __CHECKPOINT_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <complex>
#include <sys/stat.h>

#include "cpOptions.h"
#include "dataType.h"
#include "cpHelperFuncs.hpp"
#include "craftConf.h"
// ====== CP TYPES HEADERS ===== //  
#include "cpEnum.h"
#include "cpBase.hpp"
#include "cpPOD.hpp"
#include "cpArray.hpp"


#ifdef GHOST_CP 
#include "cpTypes/cpGhost/cpGhost.hpp"
#endif

#ifdef PHIST_CP 
#include "cpTypes/cpPhistMvec/cpPhistMvec.cpp"
#include "cpTypes/cpPhistSdMat/cpPhistSdMat.cpp"
#endif

#ifdef SCR
extern "C"{
	#include <scr.h>
}
#endif


#define SSTR( x ) static_cast< std::ostringstream & >( \
								        ( std::ostringstream() << std::dec << x ) ).str()


std::string exec(const char* cmd) {
  char buffer[128];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
//  if (!pipe) 
//				throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
       if (fgets(buffer, 128, pipe) != NULL)
         result = buffer;
  	}
  } 
	catch (...) {
    pclose(pipe);
  	throw;
	}
	pclose(pipe);
	result.erase(result.end()-1);
	return result;
}


class Checkpoint
{
protected:
	MPI_Comm cpMpiComm;
	std::string cpPath;
	std::string cpBasePath;	
	std::string cpPathVersion;
	std::string cpName;
	std::string cpVersionPrefix;
	bool cpUseSCR;
	bool cpCommitted;
	size_t cpVersion; 	
	size_t numBufferCps;
  bool restartStatus;

	int writeCpMetaData(const std::string filename);
//	int writeCpMetaDataSCR();
	int mkCpDir(std::string path);
	int readCpMetaData();
  int SCRreadCpMetaData();
	int deleteBackupCp();
	int PFSwrite();
	int SCRwrite();
	int PFSread();
	int SCRread();
  int getCpVersion(); 

public:
	Checkpoint(const std::string name_, const MPI_Comm cpMpiComm);
	~Checkpoint();  
	void disableSCR();
	void commit();
//  int setRestartStatus(bool toSetStatus);
//  bool getRestartStatus();
//  int readSavedCpVersion();
  bool needRestart();

	typedef std::map<const std::string,CpBase *> cp_const_map;		// TODO: check if they can be privatly declared
  typedef std::map<const std::string,CpBase *> cp_copy_map;
 	// contains points to user workspace
  	cp_const_map objects;
  	// contains asynchronous copies of objects
//  	cp_copy_map async_copies;

  int read();
 	int write();
	int update();
  
  // implementation of add() for anything that is CpBase,
  // anything else will give an error message
  int add(std::string label, CpBase * p);

  // specialized implementation of add for POD (plain old data), which is obviously
  // not derived from our base class CpBase. The same can be done for e.g. 
  // ghost_densemat
  // TODO: here a new object is created that will not be deleted unless we have some
  //       memory management in the background, the "normal" add function just adds 
  //       the pointer to the map and the person who created the object is responsible
  //       for deleting it.

  // ===== POD ===== //
  int add(std::string label, int * const i);
  int add(std::string label, float * const i);
  int add(std::string label, double * const d);
  int add(std::string label, std::complex<double> * const d);
  int add(std::string label, std::complex<float> * const d);
  // ===== POD ARRAY ===== // 
  template <class T>
  int add(std::string label, T* const arrayPtr_, const size_t nRows_);
  // ===== POD MULTI-ARRAY ===== // 
  template <class T>
  int add(std::string label, T** const arrayPtr_, const size_t nRows_, const size_t nCols_, const int toCpCol_=ALL);

  #ifdef GHOST_CP
  int add(std::string label, ghost_densemat * const GDM);
  int add(std::string label, ghost_densemat ** const GDMArray, const size_t nDenseMat_, const int toCpDenseMat_=ALL);
  //int add(std::string label, ghost_sparsemat * const GSM)
  #endif

  #ifdef PHIST_CP 
  // ===== PHIST MVEC ===== // 
  int add(std::string label, TYPE(mvec_ptr) const mvec);
  // ===== PHIST SDMAT ===== //	TODO: as MVEC, and SDMAT are both void*, they should be differenciated in some better way.  
  int add(std::string label, TYPE(sdMat_ptr) const sdMat, TYPE(sdMat_ptr) const temp);
  #endif

};

// implementation of add() for anything that is CpBase,
// anything else will give an error message
int Checkpoint::add(std::string label, CpBase * p)
{
  if( craftEnabled == 0){
		craftDbg(3, "Checkpoint::add(): CRAFT_ENABLE: %d", craftEnabled);
    return EXIT_FAILURE;
  }

  if (cpCommitted == true ){
//			printf("This checkpoint is already committed. No data can be added to checkpoint after commit() call of a checkpoint.\n");
    return EXIT_FAILURE;
  }
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
int Checkpoint::add(std::string label, int * const i){
  this->add(label, new CpPOD<int>(i, cpMpiComm));
}

int Checkpoint::add(std::string label, float * const i){
  this->add(label, new CpPOD<float>(i, cpMpiComm));
}
int Checkpoint::add(std::string label, double * const d){
  this->add(label, new CpPOD<double>(d, cpMpiComm));
}
int Checkpoint::add(std::string label, std::complex<double> * const d){
  this->add(label, new CpPOD<std::complex<double> >(d, cpMpiComm));
}	//TODO: COMPLEX write is not implemented.
int Checkpoint::add(std::string label, std::complex<float> * const d){
  this->add(label, new CpPOD<std::complex<float> >(d, cpMpiComm));
}

// ===== POD ARRAY ===== // 
template <class T>
int Checkpoint::add(std::string label, T* const arrayPtr_, const size_t nRows_){
  this->add(label, new CpArray<T>(arrayPtr_, nRows_, cpMpiComm));
}
// ===== POD MULTI-ARRAY ===== // 
template <class T>
int Checkpoint::add(std::string label, T** const arrayPtr_, const size_t nRows_, const size_t nCols_, const int toCpCol_){
  this->add(label, new CpMultiArray<T>(arrayPtr_, nRows_, nCols_, toCpCol_, cpMpiComm));
}

#ifdef GHOST_CP
// ===== GHOST DENSE MATRIX ===== //
int Checkpoint::add(std::string label, ghost_densemat * const GDM)
{
  this->add(label, new CpGhostDenseMat(GDM, cpMpiComm));
}
int Checkpoint::add(std::string label, ghost_densemat ** const GDMArray, const size_t nDenseMat_, const int toCpDenseMat_)
{
		assert (cpCommitted == false ); 
		this->add(label, new CpGhostDenseMatArray(GDMArray, nDenseMat_, toCpDenseMat_, cpMpiComm) );
}
// ===== GHOST SPARSE MATRIX ===== // TODO: add this functionality if needed by users
//int Checkpoint::add(std::string label, ghost_sparsemat * const GSM)
//{
//		this->add(label, new CpGhostSparseMat(GSM));
//}
#endif

#ifdef PHIST_CP 
// ===== PHIST MVEC ===== // 
int Checkpoint::add(std::string label, TYPE(mvec_ptr) const mvec)
{	
  this->add(label, new TYPE(CpPhistMvec)(mvec) );
}
// ===== PHIST SDMAT ===== //	TODO: as MVEC, and SDMAT are both void*, they should be differenciated in some better way.  
int Checkpoint::add(std::string label, TYPE(sdMat_ptr) const sdMat, TYPE(sdMat_ptr) const temp)
{	
  this->add(label, new TYPE(CpPhistSdMat)(sdMat, cpMpiComm) );
}
#endif

#endif
