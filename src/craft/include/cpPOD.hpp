#ifndef __CPPOD_HPP__
#define __CPPOD_HPP__

#include "cpHelperFuncs.hpp"
#include "cpBase.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

template <class T>
class CpPOD: public CpBase
{
public:
	CpPOD( T * const dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD);
	~CpPOD(){}	

	
	
private:
	T * dataPtr;
	T * asynData;
  
  int read(const std::string * filename);
	int write(const std::string * filename);
	int update();

  int readParallel(const std::string * filename);
	int readSerial(const std::string * filename);
	int writeParallel(const std::string * filename);
	int writeSerial(const std::string * filename);
	
};

template <class T>
CpPOD<T>::CpPOD( T * const dataPtr_, const MPI_Comm cpMpiComm_){
	asynData  = new T[1];
	dataPtr   = new T[1];
	dataPtr   = dataPtr_;
	*asynData = *dataPtr;
	cpMpiComm = cpMpiComm_;
//	std::cout << "ADDED POD val is: " << *dataPtr << "  " << *asynData << std::endl;
}

template <class T>
int CpPOD<T>::read(const std::string * filename)
{
#ifdef MPIIO
  readParallel(filename);
#else
  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::write(const std::string * filename)				// filename is same if MPIIO is ON else filename is different
{
#ifdef MPIIO																			// PFS checkpoints
  writeParallel(filename);
  sync(); 
#else																							// either in case of SCR or individual PFS checkpoints
  sync(); 
  writeSerial(filename);
  sync(); 
#endif
  return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::update()
{
	craftDbg(2, "CpPOD::update()");
  *asynData = *dataPtr;
  return EXIT_SUCCESS;
}


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


#endif
