#include <cstdlib>
#include <iostream>
#include <complex>
#include <iostream>

#include <mpi.h>

#include "include/cpPOD.hpp"
#include "include/cpEnum.h"
#include "include/dataType.h"

template class CpPOD<int>;
template class CpPOD<unsigned int>;
template class CpPOD<long>;
template class CpPOD<unsigned long>;
template class CpPOD<double>;
template class CpPOD<float>;
template class CpPOD<char>;
template class CpPOD<unsigned char>;
template class CpPOD<short>;
template class CpPOD<unsigned short>;
template class CpPOD<std::complex<int> >;
template class CpPOD<std::complex<float> >;
template class CpPOD<std::complex<double> >;



template <class T>
CpPOD<T>::CpPOD( T * const dataPtr_, const MPI_Comm cpMpiComm_){
	asynData  = new T[1];
	dataPtr   = new T[1];
	dataPtr   = dataPtr_;
	*asynData = *dataPtr;
	cpMpiComm = cpMpiComm_;
//		std::cout << "ADDED POD val is: " << *dataPtr << "  " << *asynData << std::endl;
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
#else																							// either in case of SCR or individual PFS checkpoints
  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::update()
{
  *asynData = *dataPtr;
//		std::cout << "AsynData is:  " << *asynData << std::endl;
  return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::readParallel(const std::string * filename){
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
	MPI_File_open(cpMpiComm, (*filename).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	os = myrank * sizeof(T)*1;			// every POD is converted to char
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_read(fh, dataPtr, 1, getMpiDataType(T), &status);
	MPI_File_close(&fh);
	printf("use--MPIIO at restart\n");
	return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::readSerial(const std::string * filename){
	printf("DONT use--MPIIO at restart\n");
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
	printf("CpPOD::writeParallel()\n");
	MPI_File fh; 
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(cpMpiComm, &myrank);
	MPI_File_open(cpMpiComm, (*filename).c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	os = myrank * sizeof(T)*1;			// every POD is converted to char
	MPI_File_seek(fh, os, MPI_SEEK_SET);
  MPI_File_write(fh, asynData, 1, getMpiDataType(T), &status);
	MPI_File_close(&fh);
//	MPI_Barrier(cpMpiComm);
	return EXIT_SUCCESS;
}

template <class T>
int CpPOD<T>::writeSerial(const std::string * filename){
	printf("CpPOD::writeSerial()\n");
	std::ofstream fstr;
	fstr.open ((*filename).c_str(), std::ios::out | std::ios::binary );	
	if(fstr.is_open()){
		fstr.write( (char *)asynData, sizeof (T) );
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
