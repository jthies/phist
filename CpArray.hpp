#ifndef __CPARRAY_H__
#define __CPARRAY_H__

#include "CpHelperFuncs.hpp"
#include "CpBase.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>

template <class T>
class CpArray: public CpBase
{
private:

	T * dataPtr;
	T * asynData;
	size_t nRows;
	void copyArray(T * const src, T * const dst, const size_t numElems);

public:
	CpArray( T * const dataPtr_, const size_t nRows_){
		nRows = nRows_;
		asynData = new T[nRows];
		dataPtr  = new T[nRows];
		dataPtr = dataPtr_;
		copyArray(dataPtr, asynData, nRows);
//	std::cout << "POD Array is added: " << std::endl;
//	print();
	}
	~CpArray(){}	

	int read(const std::string * filename)
	{
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
	
	int write(const std::string * filename)
	{
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
	
	int update()
	{
		copyArray(dataPtr, asynData, nRows);
		return EXIT_SUCCESS;
	}
	void print();
};

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


// ===== POD MULTI-ARRAY ===== // 
template <class T>
class CpMultiArray: public CpBase
{
private:

	T ** dataPtr;
	T ** asynData;
	size_t nRows;
	size_t nCols;
	int toCpCol;
	size_t cyclicCpCounter;
	void copyMultiArray(T ** const src, T ** const dst, const size_t nRows, const size_t nCols);

public:
	CpMultiArray( T ** const dataPtr_, const size_t nRows_, const size_t nCols_, const int toCpCol_=ALL){
		nRows = nRows_;
		nCols = nCols_;
		dataPtr = dataPtr_;
		asynData = new T*[nCols];
		toCpCol = toCpCol_;
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
	~CpMultiArray(){}	

	int read(const std::string * filename)
	{
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
				std::cout << "Restarting with cyclicCpCounter = " << cyclicCpCounter << std::endl;
				fstr.read( (char *)&dataPtr[cyclicCpCounter][0], sizeof (T) * nRows );
				fstr.close();
				++cyclicCpCounter;
				if( cyclicCpCounter == nCols ){
					cyclicCpCounter = 0;	
				}
			}
			else if(toCpCol>=0){
				fstr.read( (char *)&dataPtr[toCpCol][0], sizeof (T) * nRows );
				fstr.close();
			}
			else{
				std::cerr << "Not a valid toCpCol at " << __FILE__ << " " << __LINE__ << std::endl;			
				return EXIT_FAILURE;
			}
		}
		else{
			std::cerr << "Can't open file " << *filename << __FILE__ << " " << __LINE__ << std::endl;			
			return EXIT_FAILURE;
		}
		print();
		return EXIT_SUCCESS;
	}
	
	int write(const std::string * filename)
	{
		std::ofstream fstr;
		fstr.open ((*filename).c_str(), std::ios::out | std::ios::binary);	
		if(fstr.is_open()){
			if(toCpCol == ALL){
				std::cout << "toCpCol is ALL" << std::endl;
				for(int i=0; i<nCols; ++i){
					fstr.write( (char *)&(asynData[i][0]), sizeof (T) * nRows );
				}
				fstr.close();
			}
			else if(toCpCol == CYCLIC){
				std::cout << "toCpCol is cyclic" << cyclicCpCounter << std::endl;
				fstr.write( (char *)&(asynData[cyclicCpCounter][0]), sizeof (T) * nRows );
				fstr.close();
				// ===== write the metadata file for cyclicCpCounter ===== // 
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
				std::cout << "toCpCol is " << toCpCol << std::endl;
				fstr.write( (char *)&(asynData[toCpCol][0]), sizeof (T) * nRows );
				fstr.close();
			}
			else{
				std::cerr << "Not a valid toCpCol at " << __FILE__ << " " << __LINE__ << std::endl;			
				return EXIT_FAILURE;
			}
		}
		else{
			std::cerr << "Can't open file " << *filename << __FILE__ << " " << __LINE__ << std::endl;			
			return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
	}
	
	int update()
	{
		copyMultiArray(dataPtr, asynData, nRows, nCols);
		return EXIT_SUCCESS;
	}
	void print();
};

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

#endif

