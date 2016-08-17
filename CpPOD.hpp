#ifndef __CPPOD_H__
#define __CPPOD_H__

#include "CpHelperFuncs.hpp"
#include "CpBase.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

template <class T>
class CpPOD: public CpBase
{
private:

	T * dataPtr;
	T * asynData;

public:
	CpPOD( T * const dataPtr_){
		asynData = new T[1];
		dataPtr  = new T[1];
		dataPtr = dataPtr_;
		*asynData = *dataPtr;
		std::cout << "ADDED POD val is: " << *dataPtr << "  " << *asynData << std::endl;
	}
	~CpPOD(){}	

	void read(const std::string * filename)
	{
		std::ifstream fstr;
		fstr.open ((*filename).c_str(), std::ios::in | std::ios::binary);	
		fstr.read( (char*)dataPtr , sizeof(T));	
		std::cout << "dataPtr read is: " << *dataPtr << std::endl;
		fstr.close();
		return;
	}
	
	void write(const std::string * filename)
	{
		std::ofstream fstr;
		fstr.open ((*filename).c_str(), std::ios::out | std::ios::app | std::ios::binary );	
		if(fstr.is_open()){
			fstr.write( (char *)asynData, sizeof (T) );
		}
		fstr.close();
		return;
	}
	
	void update()
	{
		*asynData = *dataPtr;
		std::cout << "AsynData is:  " << *asynData << std::endl;
		return;
	}
};


#endif

