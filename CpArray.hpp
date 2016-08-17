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
		std::cout << "POD Array is added: " << std::endl;
		//print();
	}
	~CpArray(){}	

	void read(const std::string * filename)
	{
		std::ifstream fstr;
		fstr.open ((*filename).c_str());	
		std::string line;
		for(size_t i = 0; i < nRows; ++i){
			getline (fstr, line);
			asynData[i] = StringToNumber<T>(line);	// TODO: perhaps the array can be directly read into dataPtr.
			dataPtr[i] = asynData[i];
		}	
		fstr.close();
		std::cout << "POD Array after reading is: " << std::endl;
		//print();
		return;
	}
	
	void write(const std::string * filename)
	{
		std::ofstream fstr;
		fstr.open ((*filename).c_str());	
		for(size_t i = 0; i < nRows; ++i){
			fstr << std::setprecision(32) << asynData[i] << std::endl;	
		}
		fstr.close();
		return;
	}
	
	void update()
	{
		copyArray(dataPtr, asynData, nRows);
		return;
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
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << "dataPtr[" << i << "] = " << dataPtr[i] << std::endl;
	}
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << "asynData[" << i << "] = " << asynData[i] << std::endl;
	}	
	return;
}
#endif

