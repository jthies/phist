/**
 * Author: Faisal Shahzad
 * This file contains two classes for arrays and multiarrays of int, double, float datatypes.
 *
 */

#ifndef __CP_ARRAY_H__
#define __CP_ARRAY_H__

#include "enum.h"

template <class T>
class CpArray 
{
public:
		CpArray();
		~CpArray();
		int add(const std::string array_name, T * const item, const size_t rows, const MPI_Comm FtComm, const std::string cpPath_ );
		int print();
		int update();	// check it the item is defined. if it is defined, update according to nRows and array
	   	int write();	
		int read();
private:
		std::string cpPath;
		MPI_Comm cpMpiComm;
		T * array;
		T * arrayPtr;		// TODO: the pointer should be constant. But not the data, as data has to be read/written on this array after restart. it can not even be a constant pointer because the pointer has to be assigned a value in he add function.
		std::string name;
		size_t nRows;
		void copyArray(T * const src, T * const dst, const size_t numElems);
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
int CpArray<T>::add(const std::string array_name, T * const item, const size_t rows , const MPI_Comm FtComm, const std::string cpPath_)
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
	return 0;
}

template <class T>
int CpArray<T>::print(){
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << "array[" << i << "] = " << array[i] << endl;
	}
	for(size_t i = 0; i < this->nRows; ++i){
		std::cout << "arrayPtr[" << i << "] = " << arrayPtr[i] << endl;
	}	
	return 0;
}

template <class T>
int CpArray<T>::update(){
	copyArray(arrayPtr, array, nRows);
	return 0;
}

template <class T>
int CpArray<T>::write(){
	std::cout << "writing file now " << name << endl ;
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
	std::cout << "read file now" << name << endl;
	for(int i = 0; i< nRows; ++i){
		std::cout << array[i] << endl;
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
	copyArray(array, arrayPtr, nRows);	
	print();
	return 0;
}

template <class T>
void CpArray<T>::copyArray(T * const src, T * const dst, const size_t numElems){
	for(size_t i= 0; i< numElems; ++i){
		dst[i] = src[i];
	}
}	


template <class T>
class CpMulArray 
{
public:
		CpMulArray();
		~CpMulArray();
		int add(const std::string array_name, T** const item, const size_t rows , const size_t cols, const MPI_Comm FtComm, const std::string cpPath_, const int toCpCol_);
		int print();
		int update();
		int write();
		int read();

private:
		std::string cpPath;
		std::string name;
		MPI_Comm cpMpiComm;
		T ** array;
		T ** arrayPtr;
		size_t nRows;
		size_t nCols;
		int toCpCol;			// this should be read as ENUM. 0-N= specified columns, -1 = ALL, -2 = CYCLIC, starting with 0 
		void copyMulArray(T ** arrayPtr, T** array, const size_t nRows, const size_t nCols);

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
int CpMulArray<T>::add(const std::string array_name, T** const item, const size_t rows , const size_t cols, const MPI_Comm FtComm, const std::string cpPath_, const int toCpCol_)
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
	for(size_t i = 0; i < nCols; ++i){
		for(size_t j = 0; j < nRows; ++j){
				std::cout << "arrayPtr[" << i << "][" << j << "] = " << arrayPtr[i][j] << endl;
		}
	}	
	return 0;
}

template <class T>
int CpMulArray<T>::update(){ 		// TODO: 	should only update the array to be written
	copyMulArray(arrayPtr, array, nRows, nCols);	
	print();	
	return 0;
}



template <class T>
int CpMulArray<T>::write(){
	std::cout << "writing CpMulArray file now " << name << endl ;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);	
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

	FILE * fp;
	if( NULL == (fp = fopen(filename, "w+")) ) {
		fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
		return -1;
	}
	// TODO: 	write should be done according to the write-option given by the user. about which column to write
	
	if(toCpCol == ALL){
		std::cout << "writing all MULTIVEC" << endl;
		for(int i = 0; i < nCols ; ++i){
			fwrite(&array[i][0], sizeof(T), nRows, fp);
		}
	}
	if(toCpCol == CYCLIC){ 
		static size_t cyclicCpCounter = 0 ;				// this counter is responsible for determining which Col is to be checkpoinited in the cyclic case. 	
		std::cout << "writing CpCol:" << cyclicCpCounter << endl;
		fwrite(&array[cyclicCpCounter][0], sizeof(T), nRows, fp);
		cyclicCpCounter++;
		if(cyclicCpCounter == nCols){
			cyclicCpCounter = 0;
		}
	}
	if(toCpCol >= 0 ){
		std::cout << "writing toCpCol" << toCpCol<< endl;
		fwrite(&array[toCpCol][0], sizeof(T), nRows, fp);	
	}
	
	fclose(fp);

	return 0;
}

template <class T>
int CpMulArray<T>::read(){
	std::cout << "read CpMulArray file now" << name << endl;
	int myrank = -1;
	MPI_Comm_rank(cpMpiComm, &myrank);			// TODO: should be FT_comm
	char * filename = new char[256];
	sprintf(filename, "%s/%s-rank%d.cp", cpPath.c_str(), name.c_str(), myrank);

	FILE * fp;
	if( NULL == (fp = fopen(filename, "r")) ) {
		fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
		return -1;
	}
	if(toCpCol == ALL){
		std::cout << "reading all MULTIVEC" << endl;
		for(int i = 0; i < nCols ; ++i){
			fread(&array[i][0], sizeof(T), nRows, fp);
		}
	}
	if(toCpCol == CYCLIC){
		// TODO: 	check the restart with cyclic case. restart will have to determine which
		std::cerr << "Restart with CYCLIC case is not reading yet." << std::endl;
	/*	static size_t cyclicCpCounter = 0 ;				//which Col is to be read in the cyclic case. 	
		std::cout << "reading CpCol:" << cyclicCpCounter << endl;
		fread(&array[cyclicCpCounter][0], sizeof(T), nRows, fp);
		cyclicCpCounter++;
		if(cyclicCpCounter == nCols){
			cyclicCpCounter = 0;
		}
		*/
	}
	if(toCpCol >= 0 ){
		std::cout << "reading toCpCol" << toCpCol<< endl;
		fread(&array[toCpCol][0], sizeof(T), nRows, fp);	
	}	
	for(int i=0; i< nCols; ++i){
		fread( &array[i][0], sizeof(T), nRows, fp );
	}
	fclose(fp);
	copyMulArray(array, arrayPtr, nRows, nCols);	
	print();
	return 0;
}

template <class T>
void CpMulArray<T>::copyMulArray(T ** arrayPtr, T** array, const size_t nRows, const size_t nCols){
	for(size_t i = 0; i < nCols; ++i){
		for(size_t j = 0; j < nRows; ++j){
			array[i][j] = arrayPtr[i][j];
		//	std::cout << "array[" << i << "][" << j << "] = " << array[i][j] << endl;
		}
	}	
}

#endif
