#ifndef __CPMKLARRAY_HPP__
#define __CPMKLARRAY_HPP__

#include <cstdlib>
#include <fstream>
#include <iostream>
#include "cpHelperFuncs.hpp"


#include "../../cpEnum.h"
#include "../../cpBase.hpp"
#include "mkl.h"

template <class T>
static void copyMklArray(T * const src, T * const dst, const size_t nElem_){
	for(size_t i= 0; i< nElem_; ++i){
    dst[i].real = src[i].real;    
    dst[i].imag = src[i].imag;    
	}
  return;
}

// ===== MKL_Complex8 * (Array) ===== //
class CpMklComplex8Array: public CpBase  
{
public:
	CpMklComplex8Array(MKL_Complex8 *  dataPtr_, const size_t nElem_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpMklComplex8Array(){}	
private:
	MKL_Complex8 * dataPtr;
	MKL_Complex8 * asynData;
  size_t nElem; 

	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
};

CpMklComplex8Array::CpMklComplex8Array(MKL_Complex8 *  dataPtr_, const size_t nElem_, const MPI_Comm cpMpiComm_)
{
  nElem = nElem_; 
  dataPtr   = new MKL_Complex8[nElem];
  asynData  = new MKL_Complex8[nElem];
	dataPtr   = dataPtr_;
  cpMpiComm = cpMpiComm_;
  copyMklArray(dataPtr, asynData, nElem);
}

int CpMklComplex8Array::update(){
  copyMklArray(dataPtr, asynData, nElem);
  return EXIT_SUCCESS;
}

int CpMklComplex8Array::read(const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  readParallel(filename);
#else 
  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex8Array::write( const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  writeParallel(filename);
#else 
  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex8Array::writeParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex8Array::writeParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex8Array = 2*sizeof(float)*nElem;
	os = myrank * sizeofMKL_Complex8Array;
	MPI_File_seek(fh, os, MPI_SEEK_SET);
	for(size_t i= 0; i< nElem; ++i){
    MPI_File_write(fh, &(asynData[i].real), 1, getMpiDataType(float), &status);
    MPI_File_write(fh, &(asynData[i].imag), 1, getMpiDataType(float), &status);
  }
	MPI_File_close(&fh);
  sync();
	return EXIT_SUCCESS;
}

int CpMklComplex8Array::writeSerial(const std::string * filename){ // TODO: check
	craftDbg(2, "CpMklComplex8Array::writeSerial()");
	std::ofstream fstr;
   
	fstr.open ((*filename).c_str(), std::ofstream::binary );	
	if(fstr.is_open()){
	  for(size_t i= 0; i< nElem; ++i){
		  fstr.write( (char *)&(asynData[i].real), sizeof (float) );
		  fstr.write( (char *)&(asynData[i].imag), sizeof (float) );
    }
		fstr.close();
    sync();
	}
	else{
	  int myrank;
	  MPI_Comm_rank(cpMpiComm, &myrank);
		std::cerr << myrank << "Can't open file -- " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


int CpMklComplex8Array::readParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex8Array::readParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex8Array = 2*sizeof(float)*nElem; // as MKL_Complex8Array as 2 float data members.
	os = myrank * sizeofMKL_Complex8Array;	
	MPI_File_seek(fh, os, MPI_SEEK_SET);
	for(size_t i= 0; i< nElem; ++i){
    MPI_File_read(fh, &(dataPtr[i].real), 1, getMpiDataType(float), &status);
    MPI_File_read(fh, &(dataPtr[i].imag), 1, getMpiDataType(float), &status);
  }
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}


int CpMklComplex8Array::readSerial(const std::string * filename){
	craftDbg(2, "CpMklComplex8Array::readSerial()");
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
	  for(size_t i= 0; i< nElem; ++i){
		  fstr.read( (char *)&(dataPtr[i].real), sizeof (float) );
		  fstr.read( (char *)&(dataPtr[i].imag), sizeof (float) );
    }
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}



// ===== MKL_Complex16 * (Array) ===== //
class CpMklComplex16Array: public CpBase  
{
public:
	CpMklComplex16Array(MKL_Complex16 *  dataPtr_, const size_t nElem_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpMklComplex16Array(){}	
	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:
	MKL_Complex16 * dataPtr;
	MKL_Complex16 * asynData;
  size_t nElem; 

  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
};

CpMklComplex16Array::CpMklComplex16Array(MKL_Complex16 *  dataPtr_, const size_t nElem_, const MPI_Comm cpMpiComm_)
{
  nElem = nElem_; 
  dataPtr   = new MKL_Complex16[nElem];
  asynData  = new MKL_Complex16[nElem];
	dataPtr   = dataPtr_;
  cpMpiComm = cpMpiComm_;
  copyMklArray(dataPtr, asynData, nElem);
}

int CpMklComplex16Array::update(){
  copyMklArray(dataPtr, asynData, nElem);
  return EXIT_SUCCESS;
}

int CpMklComplex16Array::read(const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  readParallel(filename);
#else 
  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex16Array::write( const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  writeParallel(filename);
#else 
  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex16Array::writeParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex16Array::writeParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex8Array = 2*sizeof(double)*nElem;
	os = myrank * sizeofMKL_Complex8Array;
	MPI_File_seek(fh, os, MPI_SEEK_SET);
	for(size_t i= 0; i< nElem; ++i){
    MPI_File_write(fh, &(asynData[i].real), 1, getMpiDataType(double), &status);
    MPI_File_write(fh, &(asynData[i].imag), 1, getMpiDataType(double), &status);
  }
	MPI_File_close(&fh);
  sync();
	return EXIT_SUCCESS;
}

int CpMklComplex16Array::writeSerial(const std::string * filename){ // TODO: check
	craftDbg(2, "CpMklComplex16Array::writeSerial()");
	std::ofstream fstr;
   
	fstr.open ((*filename).c_str(), std::ofstream::binary );	
	if(fstr.is_open()){
	  for(size_t i= 0; i< nElem; ++i){
		  fstr.write( (char *)&(asynData[i].real), sizeof (double) );
		  fstr.write( (char *)&(asynData[i].imag), sizeof (double) );
    }
		fstr.close();
    sync();
	}
	else{
	  int myrank;
	  MPI_Comm_rank(cpMpiComm, &myrank);
		std::cerr << myrank << "Can't open file -- " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


int CpMklComplex16Array::readParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex16Array::readParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex8Array = 2*sizeof(double)*nElem; // as MKL_Complex8Array as 2 double data members.
	os = myrank * sizeofMKL_Complex8Array;	
	MPI_File_seek(fh, os, MPI_SEEK_SET);
	for(size_t i= 0; i< nElem; ++i){
    MPI_File_read(fh, &(dataPtr[i].real), 1, getMpiDataType(double), &status);
    MPI_File_read(fh, &(dataPtr[i].imag), 1, getMpiDataType(double), &status);
  }
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}


int CpMklComplex16Array::readSerial(const std::string * filename){
	craftDbg(2, "CpMklComplex16Array::readSerial()");
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
	  for(size_t i= 0; i< nElem; ++i){
		  fstr.read( (char *)&(dataPtr[i].real), sizeof (double) );
		  fstr.read( (char *)&(dataPtr[i].imag), sizeof (double) );
    }
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}





#endif
