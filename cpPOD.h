/**
 * Author: Faisal Shahzad
 * This file contains class for POD types. of int, 
 *
 */

#ifndef __CPPOD_H__
#define __CPPOD_H__

#include "enum.h"
#include <fstream>
#include <iostream>
#include <complex>
#include <string>
#include <typeinfo>
#include "cp.h"
#ifdef SCR
extern "C"{
	#include <scr.h>
}
#endif
using namespace std;
template <class T>
class CpPod			// TODO: check the class types. These should be child classes of the other classes. Hide them from the end user.  
{
public:
		CpPod();
		~CpPod();
		int add(const std::string podName, T * const item, const MPI_Comm FtComm, const std::string cpPath_, const bool useSCR_);
		int print();
		int update();	// check it the item is defined. if it is defined, update according to nRows and array
	   	int write();	
		int read();
private:
		std::string cpPath;
		MPI_Comm cpMpiComm;
		std::string name;
		bool useSCR;
		T pod;
		T * podPtr;		// TODO: the pointer should be constant. But not the data, as data has to be read/written on this array after restart. it can not even be a constant pointer because the pointer has to be assigned a value in he add function.
};

template <class T>
CpPod<T>::CpPod()
{
	cpMpiComm = MPI_COMM_WORLD;
	cpPath = "";
	pod = 0;
	name = "";
}

template <class T>
CpPod<T>::~CpPod()
{
}

template <class T>
int CpPod<T>::add(const std::string podName, T * const item, const MPI_Comm FtComm, const std::string cpPath_, const bool useSCR_)
{
	name 	= podName;	
	podPtr 	= item;	
	pod		= *item;
	cpMpiComm = FtComm;
	cpPath 	= cpPath_;
	useSCR 	= useSCR_;	
	return 0;
}

template <class T>
int CpPod<T>::print(){
	std::cout << "pod: " << name << " " << pod << endl;
	return 0;
}

template <class T>
int CpPod<T>::update(){
	pod = *podPtr;
	print();
	return 0;
}

template <class T>
int CpPod<T>::write(){		// TODO: check the correctness of DOUBLE COMPLEX type

	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);
	
#ifdef SCR
	if(useSCR == true){
		char * tmpFilename = new char[256];
		strcpy(tmpFilename, filename); 
		SCR_Route_file(tmpFilename, filename);
		printf("new filename %s\n", filename);
	}
#endif	
	std::ofstream fstr;// (filename);	
	fstr.open (filename);	
	fstr << name.c_str() << "\t" << pod;		
/*	if( typeid(T) == typeid(double complex) ){	// TODO: check the correctness of DOUBLE COMPLEX type
		std::cout << "ERROR: complex write is not yet tested \n" << endl;
	}*/	
	fstr.close();
	return 0;
}

template <class T>
int CpPod<T>::read(){
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

#ifdef SCR
	if(useSCR == true){
		char * tmpFilename = new char[256];
		strcpy(tmpFilename, filename); 
		SCR_Route_file(tmpFilename, filename);
		printf("new filename %s\n", filename);
	}
#endif	

	std::ifstream fstr (filename);
	char * readStr = new char[256];
	fstr >> readStr;
	name = readStr;
	
	if( typeid(T) == typeid(int) ){
		fstr >> readStr;
		pod = (T) atoi(readStr); 
		*podPtr = (T) atoi(readStr); 
	} 
	if( typeid(T) == typeid(double) ){
		fstr >> readStr;
		pod = (T) atof(readStr); 
		*podPtr = (T) atof(readStr); 
	}
	if( typeid(T) == typeid(float) ){
		fstr >> readStr;
		pod = (T) atof(readStr); 
		*podPtr = (T) atof(readStr); 
	}
/*	if( typeid(T) == typeid(complex<double> ) ){	// TODO: check the correctness of DOUBLE COMPLEX type
		fstr >> readStr;
		complex<double>  tempComplex = 0;	
		std::cout << "ERROR: complex<double>  read is not yet tested \n" << endl;
		pod = (T) atof (readStr);
		*podPtr = (T) atof (readStr);
	}*/ 
	fstr.close();
//	double tt = std::stod(tmp2);
	return 0;
}

#endif
