#ifndef __CPMKL_HPP__
#define __CPMKL_HPP__

#include <cstdlib>
#include <fstream>
#include <iostream>
#include "cpHelperFuncs.hpp"


#include "../../cpEnum.h"
#include "../../cpBase.hpp"
#include "mkl.h"

// ===== MKL_Complex8 ===== //
class CpMklComplex8: public CpBase  
{
public:
	CpMklComplex8(MKL_Complex8 *  dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpMklComplex8(){}	
	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:
	MKL_Complex8 * dataPtr;
	MKL_Complex8 * asynData;

  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
};

CpMklComplex8::CpMklComplex8(MKL_Complex8 *  dataPtr_, const MPI_Comm cpMpiComm_)
{
  dataPtr   = new MKL_Complex8[1];
  asynData  = new MKL_Complex8[1];
	dataPtr   = dataPtr_;
	*asynData = *dataPtr;
  cpMpiComm = cpMpiComm_;
}

int CpMklComplex8::update(){
  *asynData = *dataPtr;
  return EXIT_SUCCESS;
}

int CpMklComplex8::read(const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  readParallel(filename);
#else 
  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex8::write( const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  writeParallel(filename);
#else 
  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex8::writeParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex8::writeParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex8 = 2*sizeof(float);
	os = myrank * sizeofMKL_Complex8;
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_write(fh, &(asynData->real), 1, getMpiDataType(float), &status);
  MPI_File_write(fh, &(asynData->imag), 1, getMpiDataType(float), &status);
	MPI_File_close(&fh);
  sync();
	return EXIT_SUCCESS;
}

int CpMklComplex8::writeSerial(const std::string * filename){
	craftDbg(2, "CpMklComplex8::writeSerial()");
	std::ofstream fstr;
   
	fstr.open ((*filename).c_str(), std::ofstream::binary );	
	if(fstr.is_open()){
		fstr.write( (char *)&(asynData->real), sizeof (float) );
		fstr.write( (char *)&(asynData->imag), sizeof (float) );
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

int CpMklComplex8::readParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex8::readParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex8 = 2*sizeof(float); // as MKL_Complex8 as 2 float data members.
	os = myrank * sizeofMKL_Complex8;	
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_read(fh, &(dataPtr->real), 1, getMpiDataType(float), &status);
  MPI_File_read(fh, &(dataPtr->imag), 1, getMpiDataType(float), &status);
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}

int CpMklComplex8::readSerial(const std::string * filename){
	craftDbg(2, "CpMklComplex8::readSerial()");
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
		fstr.read( (char *)&(dataPtr->real), sizeof (float) );
		fstr.read( (char *)&(dataPtr->imag), sizeof (float) );
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

// ===== MKL_Complex16 ===== //
class CpMklComplex16: public CpBase  
{
public:
	CpMklComplex16(MKL_Complex16 *  dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpMklComplex16(){}	
	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:
	MKL_Complex16 * dataPtr;
	MKL_Complex16 * asynData;

  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
};

CpMklComplex16::CpMklComplex16(MKL_Complex16 *  dataPtr_, const MPI_Comm cpMpiComm_)
{
  dataPtr   = new MKL_Complex16[1];
  asynData  = new MKL_Complex16[1];
	dataPtr   = dataPtr_;
	*asynData = *dataPtr;
  cpMpiComm = cpMpiComm_;
}

int CpMklComplex16::update(){
  *asynData = *dataPtr;
  return EXIT_SUCCESS;
}

int CpMklComplex16::read(const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  readParallel(filename);
#else 
  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex16::write( const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
  writeParallel(filename);
#else 
  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklComplex16::writeParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex16::writeParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex16 = 2*sizeof(double);
	os = myrank * sizeofMKL_Complex16;
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_write(fh, &(asynData->real), 1, getMpiDataType(double), &status);
  MPI_File_write(fh, &(asynData->imag), 1, getMpiDataType(double), &status);
	MPI_File_close(&fh);
  sync();
	return EXIT_SUCCESS;
}

int CpMklComplex16::writeSerial(const std::string * filename){
	craftDbg(2, "CpMklComplex16::writeSerial()");
	std::ofstream fstr;
   
	fstr.open ((*filename).c_str(), std::ofstream::binary );	
	if(fstr.is_open()){
		fstr.write( (char *)&(asynData->real), sizeof (double) );
		fstr.write( (char *)&(asynData->imag), sizeof (double) );
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

int CpMklComplex16::readParallel(const std::string * filename){
	craftDbg(2, "CpMklComplex16::readParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
  int sizeofMKL_Complex16 = 2*sizeof(double); // as MKL_Complex16 as 2 double data members.
	os = myrank * sizeofMKL_Complex16;	
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_read(fh, &(dataPtr->real), 1, getMpiDataType(double), &status);
  MPI_File_read(fh, &(dataPtr->imag), 1, getMpiDataType(double), &status);
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}

int CpMklComplex16::readSerial(const std::string * filename){
	craftDbg(2, "CpMklComplex16::readSerial()");
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
		fstr.read( (char *)&(dataPtr->real), sizeof (double) );
		fstr.read( (char *)&(dataPtr->imag), sizeof (double) );
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


/*
// ===== MKL_INT ===== //TODO:note:not necessary as it is same as normal int or __int64, thus cpPOD class should work for this.
class CpMklInt: public CpBase  
{
public:
	CpMklInt(MKL_INT *  dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpMklInt(){}	
	int update();
	int write( const std::string * filename);
	int read(const std::string * filename);

private:
	MKL_INT * dataPtr;
	MKL_INT * asynData;

};

CpMklInt::CpMklInt(MKL_INT *  dataPtr_, const MPI_Comm cpMpiComm_){
  dataPtr   = new MKL_INT[1];
  asynData  = new MKL_INT[1];
	dataPtr   = dataPtr_;
	*asynData = *dataPtr;
  cpMpiComm = cpMpiComm_;
}

int CpMklInt::update(){
  *asynData = *dataPtr;
  return EXIT_SUCCESS;
}

int CpMklInt::read(const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
//  readParallel(filename);
#else 
//  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

int CpMklInt::write( const std::string * filename){
#ifdef MPIIO																										// Parallel PFS IO 
//  writeParallel(filename);
#else 
//  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}
*/
/*
template <class T>
int CpPOD<T>::readParallel(const std::string * filename){
	craftDbg(2, "CpPOD::readParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	os = myrank * sizeof(T)*1;			// every POD is converted to char
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_read(fh, dataPtr, 1, getMpiDataType(T), &status);
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::readSerial(const std::string * filename){
	craftDbg(2, "CpPOD::readSerial()");
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
		fstr.read( (char *) dataPtr , sizeof(T));	
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::writeParallel(const std::string * filename){
	craftDbg(2, "CpPOD::writeParallel()");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open(cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	os = myrank * sizeof(T)*1;			// every POD is converted to char
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_write(fh, asynData, 1, getMpiDataType(T), &status);
	MPI_File_close(&fh);
//	MPI_Barrier(cpMpiComm);
  sync();
	return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::writeSerial(const std::string * filename){
	craftDbg(2, "CpPOD::writeSerial()");
	std::ofstream fstr;
   
	fstr.open ((*filename).c_str(), std::ofstream::binary );	
	if(fstr.is_open()){
		fstr.write( (char *)asynData, sizeof (T) );
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
*/
#endif
