#ifndef __CPPOD_H__
#define __CPPOD_H__

#include "helperFuncs.hpp"
#include "CpBase.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>

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
		fstr.open ((*filename).c_str());	
		std::string line;
		getline (fstr, line); 
		*asynData = StringToNumber<T>(line);
		*dataPtr = *asynData;
		fstr.close();
		return;
	}
	
	void write(const std::string * filename)
	{
		std::ofstream fstr;
		fstr.open ((*filename).c_str());	
		fstr << *asynData;
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

