#ifndef __CPLIB_H__
#define __CPLIB_H__

#include <map>
#include <vector>
#include <cstdlib>
#include <string>
#include <ghost.h>
#include <ghost/types.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

//class CpArray;

template <class T>
class CpArray 
{
public:
		CpArray();
		~CpArray();
		int add(const std::string array_name, const T * const item, const size_t rows, const MPI_Comm FtComm, const std::string cpPath_ );
		int print();
		int update();	// check it the item is defined. if it is defined, update according to nRows and array
	   	int write();	
		int read();
private:
		std::string cpPath;
		MPI_Comm cpMpiComm;
		T * array;
		const T * arrayPtr;
		std::string name;
		size_t nRows;
};

template <class T>
CpArray<T>::CpArray()
{
	nRows = 0;
	cpMpiComm = MPI_COMM_WORLD;
	cpPath = "";
	array = NULL;
	name = "";
}

template <class T>
CpArray<T>::~CpArray()
{
	array = NULL;
	//arrayPtr = NULL;
}

template <class T>
int CpArray<T>::add(const std::string array_name, const T * const item, const size_t rows , const MPI_Comm FtComm, const std::string cpPath_)
{
	nRows = rows;
	arrayPtr = item;	
	name = array_name;	
	cpMpiComm = FtComm;
	cpPath = cpPath_;
	array = new T[nRows];
	for(size_t i = 0; i < nRows; ++i){
		array[i] = item[i];
	}

//	for(size_t i = 0; i < nRows; ++i){
//		cout << array[i] << endl; 
//	}
	return 0;
}

template <class T>
int CpArray<T>::print(){
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << "array[" << i << "] = " << array[i] << endl;
	}	
	return 0;
}

template <class T>
int CpArray<T>::update(){
	for(size_t i = 0; i < nRows; ++i){
		array[i] = arrayPtr[i];
	}	
	return 0;
}

template <class T>
int CpArray<T>::write(){
	cout << "writing file now " << name << endl ;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

	FILE * fp;
	if( NULL == (fp = fopen(filename, "w+")) ) {
		fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
		return -1;
	}
	fwrite(array, sizeof(T), nRows, fp);
	fclose(fp);

	return 0;
}

template <class T>
int CpArray<T>::read(){
	cout << "read file now" << name << endl;
	for(int i = 0; i< nRows; ++i){
		cout << array[i] << endl;
	}

	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

	FILE * fp;
	if( NULL == (fp = fopen(filename, "r")) ) {
		fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
		return -1;
	}
	fread( array, sizeof(T), nRows, fp );
	fclose(fp);

	return 0;
}

enum CpCol{
	ALL = -1,
	CYCLIC = -2
};

template <class T>
class CpMulArray 
{
public:
		CpMulArray();
		~CpMulArray();
		int add(const std::string array_name, const T* const* item, const size_t rows , const size_t cols, const MPI_Comm FtComm, const std::string cpPath_, const int toCpCol_);
		int print();
		int update();
		int write();
		int read();

private:
		std::string cpPath;
		std::string name;
		MPI_Comm cpMpiComm;
		T ** array;
		const T * const* arrayPtr;
		size_t nRows;
		size_t nCols;
		int toCpCol;			// this should be read as ENUM. 0-N= specified columns, -1 = ALL, -2 = CYCLIC, starting with 0 
};

template <class T>
CpMulArray<T>::CpMulArray()
{
	cpPath = "";
	name = "";
	cpMpiComm = MPI_COMM_WORLD;
	nRows = 0;
	nCols = 0;
	array = NULL;
	toCpCol = ALL;
}

template <class T>
CpMulArray<T>::~CpMulArray()
{
//	array = NULL;
}

template <class T>
int CpMulArray<T>::add(const std::string array_name, const T* const* item, const size_t rows , const size_t cols, const MPI_Comm FtComm, const std::string cpPath_, const int toCpCol_)
{
	cpPath = cpPath_;
	name = array_name;
	nRows = rows;	
	nCols = cols;
	cpMpiComm = FtComm;
	arrayPtr = item;	

	array = new T*[nCols];
	for(size_t i = 0; i < nCols; ++i){
		array[i] = new T[nRows];
	}
	for(size_t i = 0; i < nCols; ++i){
		for(size_t j = 0; j < nRows; ++j){
			array[i][j] = item[i][j];
		}
	}
	toCpCol = toCpCol_;
	return 0;
}

template <class T>
int CpMulArray<T>::print(){
	for(size_t i = 0; i < nCols; ++i){
		for(size_t j = 0; j < nRows; ++j){
				std::cout << "array[" << i << "][" << j << "] = " << array[i][j] << endl;
		}
	}	
	return 0;
}

template <class T>
int CpMulArray<T>::update(){
	for(size_t i = 0; i < nCols; ++i){
		for(size_t j = 0; j < nRows; ++j){
			array[i][j] = arrayPtr[i][j];
			std::cout << "array[" << i << "][" << j << "] = " << array[i][j] << endl;
		}
	}	
	return 0;
}



template <class T>
int CpMulArray<T>::write(){
	cout << "writing CpMulArray file now " << name << endl ;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

	FILE * fp;
	if( NULL == (fp = fopen(filename, "w+")) ) {
		fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
		return -1;
	}
	// TODO: 	write should be done according to the write-option given by the user. about which column to write
	
	if(toCpCol == ALL){
		cout << "writing all MULTIVEC" << endl;
		for(int i = 0; i < nCols ; ++i){
			fwrite(&array[i][0], sizeof(T), nRows, fp);
		}
	}
	if(toCpCol == CYCLIC){ 
		static size_t CpCol = 0 ;				// this counter is responsible for determining which Col is to be checkpoinited in the cyclic case. 	
		cout << "writing CpCol:" << CpCol << endl;
		fwrite(&array[CpCol][0], sizeof(T), nRows, fp);
		CpCol ++;
		if(CpCol == nCols){
			CpCol = 0;
		}
	}
	if(toCpCol >= 0 ){
		cout << "writing toCpCol" << toCpCol<< endl;
		fwrite(&array[toCpCol][0], sizeof(T), nRows, fp);	
	}
	
	fclose(fp);

	return 0;
}

template <class T>
int CpMulArray<T>::read(){
	cout << "read CpMulArray file now" << name << endl;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

	FILE * fp;
	if( NULL == (fp = fopen(filename, "r")) ) {
		fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
		return -1;
	}
	// TODO: 	read should be done according to the write-option given by the user. about which column to write
	for(int i=0; i< nCols; ++i){
		fread( &array[i][0], sizeof(T), nRows, fp );
	}
	fclose(fp);

	return 0;
}





class CP{
private:
	MPI_Comm cpmpicomm;
	std::string cpPath;
	bool cpCommitted;
	// ===== POD ===== 
	std::map<const std::string, const int * const> intPod;
	std::map<const std::string, int > intPodAsync;

	std::map<const std::string, const double * const> doublePod;
	std::map<const std::string, double > doublePodAsync;

	std::map<const std::string, const float * const > floatPod;
	std::map<const std::string, float > floatPodAsync;

	// ===== GHOST-VECTORS ===== 
	std::map<std::string, ghost_densemat *> vec;
	std::map<std::string, ghost_densemat *> vecAsync;

	ghost_densemat *vector_async;

	// ===== ARRAYS ===== 

	std::map<const std::string, CpArray<int> * > CpIntArrayMap;
	std::map<const std::string, CpArray<double> * > CpDoubleArrayMap;
	std::map<const std::string, CpMulArray<int> * > CpIntMulArrayMap;
	std::map<const std::string, CpMulArray<double> * > CpDoubleMulArrayMap;

	int CP_ADD_POD_INT		(const std::string key, const int * const value);
	int CP_ADD_POD_DOUBLE	(const std::string key, const double * const value); 
	int CP_ADD_POD_FLOAT	(const std::string key, const float * const value); 

	int CP_ADD_ARRAY_INT(const std::string key, const int * const val_array, const size_t array_size );
	int CP_ADD_ARRAY_DOUBLE(const std::string key,  const double * const val_array, const size_t array_size );

	int CP_ADD_MULTIARRAY_INT	(const std::string key, const int* const* ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );
	int CP_ADD_MULTIARRAY_DOUBLE(const std::string key, const double* const* ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );

public:
	CP();
	~CP();
	int setCpPath(const std::string cpPath_);
	int setComm(const MPI_Comm FT_Comm);
	void commit();

	int CP_ADD_VEC(const std::string key, ghost_densemat * value);
	int CP_ADD_MULTI_VEC(const std::string key, ghost_densemat ** value, int num_vecs);
	
	int CP_ADD_POD(const std::string key, const int * const val_ptr);
	int CP_ADD_POD(const std::string key, const double * const val_ptr);
	int CP_ADD_POD(const std::string key, const float * const val_ptr);

	int CP_ADD_ARRAY(const std::string key,  const int 	* const val_array, const size_t array_size);
	int CP_ADD_ARRAY(const std::string key,  const double * const val_array, const size_t array_size);

	int CP_ADD_MULTIARRAY(const std::string key, const int* const* ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);
	int CP_ADD_MULTIARRAY(const std::string key, const double* const* ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);

	int updateCp();	// this function should update the values of CP
	int writeCp();		// this function should write the updated CP
	int readCp();

	int print();
};





/*
template <class T>
class mypair {
	T a, b;
public:
	mypair (T first, T second)
    {a=first; b=second;}
   T getmax ();
};

*/

#endif



// int * getValueOfKey(const std::string &key) const ;

/*
template <class T>
T mypair<T>::getmax ()
{
	T retval;
	retval = a>b? a : b;
	return retval;
}
*/


