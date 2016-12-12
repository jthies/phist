#include <cstdlib>
#include <iostream>

#include <mpi.h>

#include "include/cpArray.hpp"
#include "include/cpEnum.h"
#include "include/dataType.h"

template class CpArray<int>;
template class CpArray<unsigned int>;
template class CpArray<long>;
template class CpArray<unsigned long>;
template class CpArray<double>;
template class CpArray<float>;
template class CpArray<char>;
template class CpArray<unsigned char>;
template class CpArray<short>;
template class CpArray<unsigned short>;

template class CpMultiArray<int>;
template class CpMultiArray<unsigned int>;
template class CpMultiArray<long>;
template class CpMultiArray<unsigned long>;
template class CpMultiArray<double>;
template class CpMultiArray<float>;
template class CpMultiArray<char>;
template class CpMultiArray<unsigned char>;
template class CpMultiArray<short>;
template class CpMultiArray<unsigned short>;



template <class T>
CpArray<T>::CpArray( T * const dataPtr_, const size_t nRows_, const MPI_Comm cpMpiComm_)
{
  nRows = nRows_;
  asynData = new T[nRows];
  dataPtr  = new T[nRows];
  dataPtr = dataPtr_;
  cpMpiComm = cpMpiComm_;
  copyArray(dataPtr, asynData, nRows);
//	std::cout << "POD Array is added: " << std::endl;
//	print();
}


template <class T>
int CpArray<T>::read(const std::string * filename)
{
#ifdef MPIIO
  readParallel(filename);
#else
  readSerial(filename);
#endif
  return EXIT_SUCCESS;
}

template <class T>
int CpArray<T>::write(const std::string * filename)
{
#ifdef MPIIO
  writeParallel(filename);
#else
  writeSerial(filename);
#endif
  return EXIT_SUCCESS;
}

template <class T>
int CpArray<T>::update()
{
  copyArray(dataPtr, asynData, nRows);
  return EXIT_SUCCESS;
}

template <class T>
int CpArray<T>::readSerial(const std::string * filename){
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
		fstr.read( (char*)dataPtr , sizeof(T) * nRows);	
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
//	std::cout << "POD Array after reading is: " << std::endl;
//	print();
	return EXIT_SUCCESS;
}

template <class T>
int CpArray<T>::readParallel(const std::string * filename){
	craftDbg(2, "CpArray<T>::readParallel");
	MPI_File fh;
	MPI_Status status;
	int myrank=-1; 
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open( cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	os = myrank * sizeof (T) * nRows;
	MPI_File_seek( fh , os, MPI_SEEK_SET);
	MPI_File_read( fh, dataPtr, nRows, getMpiDataType(T), &status); 
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}

template <class T>
int CpArray<T>::writeSerial(const std::string * filename){
	std::ofstream fstr;
	fstr.open ((*filename).c_str(), std::ios::out | std::ios::binary);	
	if(fstr.is_open()){
		fstr.write( (char *)asynData, sizeof (T) * nRows );
		fstr.close();
	}
	else{
		std::cerr << "Can't open file " << *filename << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}	

template <class T>
int CpArray<T>::writeParallel(const std::string * filename){
	craftDbg(2, "CpArray<T>::writeParallel");
	MPI_File fh;
	MPI_Status status;
	int myrank=-1; 
	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open( cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	os = myrank * sizeof (T) * nRows;
	MPI_File_seek( fh , os, MPI_SEEK_SET);
	MPI_File_write( fh, asynData, nRows, getMpiDataType(T), &status); 
	MPI_File_close(&fh);
	return EXIT_SUCCESS;
}

template <class T>
void CpArray<T>::copyArray(T * const src, T * const dst, const size_t numElems){
	for(size_t i= 0; i< numElems; ++i){
		dst[i] = src[i];
	}
}

template <class T>
void CpArray<T>::print(){
	std::cout << "dataPtr:" << std::endl;
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << dataPtr[i] << std::endl;
	}
	std::cout << "asynData:" << std::endl;
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << asynData[i] << std::endl;
	}	
	return;
}


template <class T>
CpMultiArray<T>::CpMultiArray( T ** const dataPtr_, const size_t nRows_, const size_t nCols_, const int toCpCol_, const MPI_Comm cpMpiComm_){
		nRows = nRows_;
		nCols = nCols_;
		dataPtr = dataPtr_;
		asynData = new T*[nCols];
		toCpCol = toCpCol_;
		cpMpiComm = cpMpiComm_;
		cyclicCpCounter = 0;
		for(size_t i = 0; i < nCols; ++i){
			asynData[i] = new T[nRows];
		}
		for(size_t i = 0; i < nCols; ++i){
			for(size_t j = 0; j < nRows; ++j){
				asynData[i][j] = dataPtr[i][j];
			}
		}
		copyMultiArray(dataPtr, asynData, nRows, nCols);
//	std::cout << "POD Array is added: " << std::endl;
//	print();
	}

template <class T>
int CpMultiArray<T>::read(const std::string * filename)
{
#ifdef MPIIO
	readParallel(filename);
#else
	readSerial(filename);
#endif
	return EXIT_SUCCESS;
}
	
template <class T>
int CpMultiArray<T>::write(const std::string * filename)
{
#ifdef MPIIO
	writeParallel(filename);
#else
	writeSerial(filename);
#endif
	return EXIT_SUCCESS;
}

template <class T>
int CpMultiArray<T>::update()
{
	copyMultiArray(dataPtr, asynData, nRows, nCols);
	return EXIT_SUCCESS;
}

template <class T>
int CpMultiArray<T>::readParallel(const std::string * filename)
{
	craftDbg(2, "CpMultiArray<T>::readParallel\n");
	MPI_File fh;
	MPI_Status status;
	int myrank=-1; 

	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open( cpMpiComm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	{
		if(toCpCol == ALL){
			craftDbg(3, "read toCpCol is ALL");
			os = myrank * sizeof (T) * nRows * nCols;
			MPI_File_seek( fh , os, MPI_SEEK_SET);
			for(int i = 0; i< nCols ; ++i){
				MPI_File_read( fh, &(dataPtr[i][0]), nRows, getMpiDataType(T), &status); 
			}
			MPI_File_close(&fh);
		}
		else if(toCpCol == CYCLIC){
			std::string filenameMD;
			filenameMD = *filename + ".metadata";
			std::ifstream fstrMD;
			fstrMD.open ((filenameMD).c_str(), std::ios::in | std::ios::binary );	
			if(fstrMD.is_open()){
				fstrMD.read( (char *)&(cyclicCpCounter), sizeof (size_t) );
				fstrMD.close();
			}
			craftDbg(3, "Restarting with cyclicCpCounter = %d", cyclicCpCounter);

			os = myrank * sizeof (T) * nRows;
			MPI_File_seek( fh , os, MPI_SEEK_SET);
			MPI_File_read( fh, &(dataPtr[cyclicCpCounter][0]), nRows, getMpiDataType(T), &status); 
			MPI_File_close(&fh); 
			++cyclicCpCounter;
			if( cyclicCpCounter == nCols ){
				cyclicCpCounter = 0;	
			}
		}
		else if(toCpCol>=0){
			craftDbg(3, "read toCpCol is %d", toCpCol);
			os = myrank * sizeof (T) * nRows;
			MPI_File_seek( fh , os, MPI_SEEK_SET);
			MPI_File_read( fh, &(dataPtr[toCpCol][0]), nRows, getMpiDataType(T), &status); 
			MPI_File_close(&fh);
		}
		else{
			std::cerr << "not a valid toCpCol at " << __FILE__ << " " << __LINE__ << std::endl;			
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}

template <class T>
int CpMultiArray<T>::readSerial(const std::string * filename){
	std::ifstream fstr;
	fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
	if(fstr.is_open()){
		if(toCpCol == ALL){
			for(int i=0; i<nCols; ++i){
				fstr.read( (char *)&dataPtr[i][0], sizeof (T) * nRows );
			}
			fstr.close();
		}
		else if(toCpCol == CYCLIC){ 			//TODO: testing not performed
			// ===== read the metadata file for cyclicCpCounter ===== // 
			std::string filenameMD;
			filenameMD = *filename + ".metadata";
			std::ifstream fstrMD;
			fstrMD.open ((filenameMD).c_str(), std::ios::in | std::ios::binary );	
			if(fstrMD.is_open()){
				fstrMD.read( (char *)&(cyclicCpCounter), sizeof (size_t) );
				fstrMD.close();
			}
			craftDbg(3, "Restarting with cyclicCpCounter = %d", cyclicCpCounter);
			fstr.read( (char *)&dataPtr[cyclicCpCounter][0], sizeof (T) * nRows );
			fstr.close();
			++cyclicCpCounter;
			if( cyclicCpCounter == nCols ){
				cyclicCpCounter = 0;	
			}
		}
		else if(toCpCol >=0){
			fstr.read( (char *)&dataPtr[toCpCol][0], sizeof (T) * nRows );
			fstr.close();
		}
		else{
			std::cerr << "not a valid tocpcol at " << __FILE__ << " " << __LINE__ << std::endl;			
			return EXIT_FAILURE;
		}
	}
	else{
		std::cerr << "can't open file " << *filename << __FILE__ << " " << __LINE__ << std::endl;			
		return EXIT_FAILURE;
	}
	print();
	return EXIT_SUCCESS;
}

template <class T>
int CpMultiArray<T>::writeParallel(const std::string * filename){
	craftDbg(2, "CpMultiArray<T>::writeParallel");
	MPI_File fh;
	MPI_Status status;
	int myrank=-1; 

	MPI_Comm_rank(cpMpiComm, &myrank);
  char * fname = new char[256];
  sprintf(fname, "%s",  (*filename).c_str());
	MPI_File_open( cpMpiComm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_Offset os;
	{
		if(toCpCol == ALL){
			craftDbg(3, "toCpCol is ALL");
			os = myrank * sizeof (T) * nRows * nCols;
			MPI_File_seek( fh , os, MPI_SEEK_SET);
			for(int i = 0; i< nCols ; ++i){
				MPI_File_write( fh, &(asynData[i][0]), nRows, getMpiDataType(T), &status); 
			}
			MPI_File_close(&fh);
		}
		else if(toCpCol == CYCLIC){
			craftDbg(3, "toCpCol is CYCLIC: %d", cyclicCpCounter);
			os = myrank * sizeof (T) * nRows;
			MPI_File_seek( fh , os, MPI_SEEK_SET);
			MPI_File_write( fh, &(asynData[cyclicCpCounter][0]), nRows, getMpiDataType(T), &status); 
			MPI_File_close(&fh);
			// ===== write the metadata file for cycliccpcounter ===== // 
			if(myrank==0){
				std::string filenameMD;
				filenameMD = *filename + ".metadata";
				std::ofstream fstrMD;
				fstrMD.open ((filenameMD).c_str(), std::ios::out | std::ios::binary );	
				if(fstrMD.is_open()){
					fstrMD.write( (char *)&(cyclicCpCounter), sizeof (size_t) );
					fstrMD.close();
				}
			}
			MPI_Barrier(cpMpiComm);
			++cyclicCpCounter;
			if( cyclicCpCounter == nCols ){
				cyclicCpCounter = 0;	
			}
		}
		else if(toCpCol>=0){
			craftDbg(3, "toCpCol is: %d", toCpCol);
			os = myrank * sizeof (T) * nRows;
			MPI_File_seek( fh , os, MPI_SEEK_SET);
			MPI_File_write( fh, &(asynData[toCpCol][0]), nRows, getMpiDataType(T), &status); 
			MPI_File_close(&fh);
		}
		else{
			std::cerr << "not a valid toCpCol at " << __FILE__ << " " << __LINE__ << std::endl;			
			return EXIT_FAILURE;
		}
	}

		return EXIT_SUCCESS;
}

template <class T>
int CpMultiArray<T>::writeSerial(const std::string * filename){
	std::ofstream fstr;
	fstr.open ((*filename).c_str(), std::ios::out | std::ios::binary);	
	if(fstr.is_open()){
		if(toCpCol == ALL){
		  craftDbg(3, "toCpCol is ALL");
			for(int i=0; i<nCols; ++i){
				fstr.write( (char *)&(asynData[i][0]), sizeof (T) * nRows );
			}
			fstr.close();
		}
		else if(toCpCol == CYCLIC){
			craftDbg(3, "toCpCol is CYCLIC: %d", cyclicCpCounter);
			fstr.write( (char *)&(asynData[cyclicCpCounter][0]), sizeof (T) * nRows );
			fstr.close();
			// ===== write the metadata file for cycliccpcounter ===== // 
			std::string filenameMD;
			filenameMD = *filename + ".metadata";
			std::ofstream fstrMD;
			fstrMD.open ((filenameMD).c_str(), std::ios::out | std::ios::binary );	
			if(fstrMD.is_open()){
				fstrMD.write( (char *)&(cyclicCpCounter), sizeof (size_t) );
				fstrMD.close();
			}
			++cyclicCpCounter;
			if( cyclicCpCounter == nCols ){
				cyclicCpCounter = 0;	
			}
		}
		else if(toCpCol>=0){
			craftDbg(3, "toCpCol is: %d", toCpCol);
			fstr.write( (char *)&(asynData[toCpCol][0]), sizeof (T) * nRows );
			fstr.close();
		}
		else{
			std::cerr << "not a valid toCpCol at " << __FILE__ << " " << __LINE__ << std::endl;			
			return EXIT_FAILURE;
		}
	}
	else{
		std::cerr << "can't open file " << *filename << __FILE__ << " " << __LINE__ << std::endl;			
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

template <class T>
void CpMultiArray<T>::print(){
	std::cout << "dataPtr:" << std::endl;
	for(size_t i = 0; i < this->nRows; ++i){
		for(size_t j = 0; j < this->nCols; ++j){
			std::cout << dataPtr[j][i] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << "asynData:" << std::endl;
	for(size_t i = 0; i < this->nRows; ++i){
		for(size_t j = 0; j < this->nCols; ++j){
			std::cout << asynData[j][i] << "\t";
		}
		std::cout << std::endl;
	}
	return;
}

template <class T>
void CpMultiArray<T>::copyMultiArray(T ** src, T** dst, const size_t nRows, const size_t nCols){
	for(size_t i = 0; i < nCols; ++i){
		for(size_t j = 0; j < nRows; ++j){
			dst[i][j] = src[i][j];
		}
	}	
}

