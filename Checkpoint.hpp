#ifndef __CHECKPOINT_HPP__
#define __CHECKPOINT_HPP__

#include "cp_options.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex>

// ====== CP TYPES HEADERS ===== //  
#include "CpEnum.h"
#include "CpBase.hpp"
#include "CpPOD.hpp"
#include "CpArray.hpp"
#include <sys/stat.h>

#ifdef GHOST_CP 
#include "CpGhost.hpp"
#endif

#ifdef PHIST_CP 
#include "cpTypes/cpPhistMvec/CpPhistMvec.cpp"
#include "cpTypes/cpPhistSdMat/CpPhistSdMat.cpp"
#endif

#ifdef SCR
extern "C"{
	#include <scr.h>
}
#endif


#define SSTR( x ) static_cast< std::ostringstream & >( \
								        ( std::ostringstream() << std::dec << x ) ).str()


class Checkpoint
{
protected:
	MPI_Comm cpMpiComm;
	std::string cpPath;
	std::string cpBasePath;	
	std::string cpPathVersion;
	std::string name;
	std::string cpVersionPrefix;
	bool useSCR;
	bool cpCommitted;
	size_t cpVersion; 	
	size_t numBufferCps;
	int writeCpMetaData();
	int mkCpDir(std::string path);
	int readCpMetaData();
	int deleteBackupCp();
public:
	Checkpoint(const std::string name_, const std::string cpBasePath_, const MPI_Comm cpMpiComm_);
	~Checkpoint();  
	void setCpPath (const std::string cpPath_);
	void setComm (const MPI_Comm cpMpiComm_);
	void disableSCR();
	void commit();

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
void add(std::string label, CpBase * p)
{
		assert (cpCommitted == false );
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
void add(std::string label, int * const i){
		this->add(label, new CpPOD<int>(i));
}

void add(std::string label, float * const i){
		this->add(label, new CpPOD<float>(i));
}
void add(std::string label, double * const d){
		this->add(label, new CpPOD<double>(d));
}
void add(std::string label, std::complex<double> * const d){
		this->add(label, new CpPOD<std::complex<double> >(d));
}	//TODO: COMPLEX write is not implemented.
void add(std::string label, std::complex<float> * const d){
		this->add(label, new CpPOD<std::complex<float> >(d));
}

// ===== POD ARRAY ===== // 
template <class T>
void add(std::string label, T* const arrayPtr_, const size_t nRows_){
		this->add(label, new CpArray<T>(arrayPtr_, nRows_));
}
// ===== POD MULTI-ARRAY ===== // 
template <class T>
void add(std::string label, T** const arrayPtr_, const size_t nRows_, const size_t nCols_, const int toCpCol_=ALL){
		this->add(label, new CpMultiArray<T>(arrayPtr_, nRows_, nCols_, toCpCol_));
}

#ifdef GHOST_CP
// ===== GHOST DENSE MATRIX ===== //
void add(std::string label, ghost_densemat * const GDM)
{
		this->add(label, new CpGhostDenseMat(GDM));
}
//void add(std::string label, ghost_densemat ** const GDMArray, const size_t nDenseMat_, const int toCpDenseMat_=ALL)
//{
//		assert (cpCommitted == false ); 
//		this->add(label, new CpGhostDenseMatArray(GDMArray, nDenseMat_, toCpDenseMat_) );
//}
// ===== GHOST SPARSE MATRIX ===== // TODO: add this functionality if needed by users
//void add(std::string label, ghost_sparsemat * const GSM)
//{
//		this->add(label, new CpGhostSparseMat(GSM));
//}
#endif



#ifdef PHIST_CP 
// ===== PHIST MVEC ===== // 
void add(std::string label, TYPE(mvec_ptr) const Mvec)
{	
		this->add(label, new TYPE(CpPhistMvec)(Mvec) );
}
// ===== PHIST SDMAT ===== //	TODO: as MVEC, and SDMAT are both void*, they should be differenciated in some better way.  
void add(std::string label, TYPE(sdMat_ptr) const sdMat, TYPE(sdMat_ptr) const temp)
{	
		this->add(label, new TYPE(CpPhistSdMat)(sdMat) );
}
#endif

};

Checkpoint::Checkpoint(const std::string name_, const std::string cpBasePath_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD){
	name 				= name_;
	cpMpiComm 	= cpMpiComm_;	
	cpBasePath 	= cpBasePath_;
	cpPath 			= cpBasePath_ + "/" + name;
	cpVersionPrefix = "CP-";
	cpCommitted = false;
	cpVersion 	= 0;
	numBufferCps=2;
	std::cout << "constructor is being called here" << std::endl;
	mkCpDir(cpPath);	
#ifdef SCR																		// check if CPAFTLIB was compiled with SCR
	useSCR 		= true;
#else
	useSCR 		= false;
#endif
}	

Checkpoint::~Checkpoint(){
#ifdef SCR
	if(useSCR == true){
		printf("====== SCR_Finalize is done ====== \n");
		SCR_Finalize();
	}
#endif
}

void Checkpoint::disableSCR(){
#ifdef SCR
	printf("SCR is disabled. Use SCR to reduce the impact of Checkpoint/restart\n");
	useSCR = false;
#endif
}

void Checkpoint::commit(){
	cpCommitted = true;	
}

int Checkpoint::mkCpDir(std::string path){
	std::string cmd = "mkdir -p " + path;
	system (cmd.c_str());
	return 0;
}

int Checkpoint::writeCpMetaData(){			// writes a tmp file first and then replaces it after complete write
	int myrank = -1;	
	MPI_Comm_rank(cpMpiComm, &myrank);
	if(myrank == 0){
		std::string filename, filenameTmp;
		filenameTmp = cpPath + "/" + "metadata.ckpt.tmp"; 
		filename = cpPath + "/" + "metadata.ckpt"; 
		std::ofstream fstr;
		fstr.open ((filenameTmp).c_str(), std::ios::out );	
		if(fstr.is_open()){
			fstr << cpVersion << std::endl;
			fstr.close();
		}
		std::string cmd = "mv " + filenameTmp + " " + filename;
		system (cmd.c_str());
	}
	return 0;
}

int Checkpoint::readCpMetaData(){
	std::string filename, line;
	filename = cpPath + "/" + "metadata.ckpt"; 
	std::ifstream fstr;
	fstr.open ((filename).c_str(), std::ios::in );	
	if(fstr.is_open()){
		while ( getline (fstr,line) )
		{
			cpVersion = atoi(line.c_str());
		}
		fstr.close();
		std::cout << "cpVersion: " << cpVersion << std::endl;
	}
	return 0;
}
int Checkpoint::deleteBackupCp(){
	int myrank = -1;	
	MPI_Comm_rank(cpMpiComm, &myrank);
	if(myrank == 0){
		std::string toRmDir = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion-numBufferCps-1);
		std::string cmd = "rm -r " + toRmDir;	
		struct stat sb;
		if (stat(toRmDir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
		{
			printf("EXists%d \n", cpVersion - numBufferCps - 1 );
			system ( cmd.c_str());
 		}
			else
		{
			printf("doNOTExists%d\n", cpVersion - numBufferCps - 1 );
		}
	}
	return 0;
}

#endif
