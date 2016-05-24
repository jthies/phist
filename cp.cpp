#include "cp.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <stdlib.h>
#include <cstdlib>

#define SAFE_INSERT(a) \
	if (a.second==false){ \
		printf("MAP member exists at %s (%d) \n", __FILE__, __LINE__);\
		exit(0);\
	}

CP::CP(){
	cpmpicomm = MPI_COMM_WORLD;	
	cpPath = "";
	cpCommitted = false;
}

CP::~CP(){

}

int CP::setCpPath(const std::string cpPath_){
	cpPath = cpPath_;
	return 0;
}
	
int CP::setComm(MPI_Comm FT_Comm){
	cpmpicomm = FT_Comm;
	return 0;
}

void CP::commit(){
	cpCommitted = true;	
}


// ========== GHOST-VECTOR CALLS ========== //
int CP::CP_ADD_GHOST_DENSEMAT(const std::string key, ghost_densemat * const value){	
	assert (cpCommitted == false);
	CpGhostDenseMat * cpdata = new CpGhostDenseMat[1];
	cpdata->add(key, value , cpmpicomm, cpPath);
	cpGhostDenseMat.insert(std::pair<const std::string, CpGhostDenseMat * >(key, cpdata));	
	return 0;
}


int CP::CP_ADD_GHOST_DENSEMAT_ARRAY(const std::string key, ghost_densemat ** const value, const int nVecs, const int toCpVec_)
{
	assert (cpCommitted == false);
	CpGhostDenseMatArray * cpdata= new CpGhostDenseMatArray[1];
	cpdata->add(key, value , nVecs, cpmpicomm, cpPath, toCpVec_);
	cpGhostDenseMatArray.insert( std::pair< const std::string, CpGhostDenseMatArray * >(key, cpdata));
	return 0;
}
	

// ========== POD CALLS ========== //

int CP::CP_ADD_POD(const std::string key, const int * const val_ptr){
	assert (cpCommitted == false);
	const	int * const a = val_ptr;
	 //a = val_ptr;
	CP_ADD_POD_INT(key, val_ptr);
	return 0;
}

int CP::CP_ADD_POD(std::string key, const double * const val_ptr){
	assert (cpCommitted == false);
	CP_ADD_POD_DOUBLE(key, val_ptr);	
	return 0;
}

int CP::CP_ADD_POD(std::string key, const float * const val_ptr){
	assert (cpCommitted == false)	;
	CP_ADD_POD_FLOAT(key, val_ptr);	
	return 0;
}

int CP::CP_ADD_POD_INT(const std::string key, const int * const value){
	int  a = *value;	
	SAFE_INSERT( intPod.insert(std::pair<std::string, const int * const>(key, value)) );
	SAFE_INSERT( intPodAsync.insert( std::pair<std::string, int > (key, *value)) );
	return 0;
}

int CP::CP_ADD_POD_DOUBLE(const std::string key, const double * const value){	 
	SAFE_INSERT( doublePod.insert( std::pair<std::string, const double * const >(key, value)) );
	SAFE_INSERT( doublePodAsync.insert(std::pair<std::string, double > (key, *value)) );
	return 0;
}

int CP::CP_ADD_POD_FLOAT(const std::string key, const float * const value){ 
	SAFE_INSERT( floatPod.insert( std::pair<std::string, const float * const >(key, value)) );
	SAFE_INSERT( floatPodAsync.insert( std::pair<std::string, float > (key, *value)) );
	return 0;
}


// ========== SINGLE ARRAY CALLS ========== //

int CP::CP_ADD_ARRAY(const std::string key, const int * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_INT(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY(const std::string key, const double * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_DOUBLE(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY(const std::string key, const float * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_FLOAT(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY_INT(const std::string key, const int * const val_array, const size_t array_size ){
	printf("Type of array is int\n");

	CpArray<int> * arraydata = new CpArray<int>[1];
	arraydata->add(key, val_array, array_size, cpmpicomm, cpPath );		arraydata->print();	
	arraydata->print();	
	cpIntArrayMap.insert(std::pair<const std::string, CpArray<int> * >(key, arraydata));	

	return 0;	
}

int CP::CP_ADD_ARRAY_DOUBLE(const std::string key, const double * const val_array, const size_t array_size ){
	printf("Type of array is double\n");
	CpArray<double> * arraydata = new CpArray<double>[1];
	arraydata->add(key, val_array, array_size , cpmpicomm, cpPath);
//	arraydata->print();
	cpDoubleArrayMap.insert(std::pair<const std::string, CpArray<double> * >(key, arraydata));	

	return 0;	
}

int CP::CP_ADD_ARRAY_FLOAT(const std::string key, const float * const val_array, const size_t array_size ){
	printf("Type of array is float\n");
	CpArray<float> * arraydata = new CpArray<float>[1];
	arraydata->add(key, val_array, array_size , cpmpicomm, cpPath);
//	arraydata->print();
	cpFloatArrayMap.insert(std::pair<const std::string, CpArray<float> * >(key, arraydata));	

	return 0;	
}

// ========== MULTI ARRAY CALLS ========== //

int CP::CP_ADD_MULTIARRAY(const std::string key, const int * const* ptr, const size_t nRows, const size_t nCols, const int toCpCol_ ){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_INT( key,  ptr, nRows, nCols, toCpCol_);
	return 0;
}

int CP::CP_ADD_MULTIARRAY(const std::string key, const double* const* ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_DOUBLE( key,  ptr, nRows, nCols, toCpCol_);
	return 0;
}

int CP::CP_ADD_MULTIARRAY(const std::string key, const float* const* ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_FLOAT( key,  ptr, nRows, nCols, toCpCol_);
	return 0;
}

int CP::CP_ADD_MULTIARRAY_INT(const std::string key, const int* const* ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	
	CpMulArray<int> * arraydata = new CpMulArray<int>[1];
	arraydata->add(key, ptr, nRows, nCols, cpmpicomm, cpPath, toCpCol_);	
	cpIntMulArrayMap.insert(std::pair<const std::string, CpMulArray<int> * > (key, arraydata));

	arraydata->print();
	return 0;
}

int CP::CP_ADD_MULTIARRAY_DOUBLE(const std::string key, const double * const* ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	printf ("Adding DOUBLE MULTIARRAY\n");	

	CpMulArray<double> * arraydata = new CpMulArray<double>[1];
	arraydata->add(key, ptr, nRows, nCols, cpmpicomm, cpPath, toCpCol_);	
	cpDoubleMulArrayMap.insert(std::pair<const std::string, CpMulArray<double> * > (key, arraydata));
	
	arraydata->print();	
	return 0;
}


int CP::CP_ADD_MULTIARRAY_FLOAT(const std::string key, const float * const* ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	printf ("Adding FLOAT MULTIARRAY\n");	

	CpMulArray<float> * arraydata = new CpMulArray<float>[1];
	arraydata->add(key, ptr, nRows, nCols, cpmpicomm, cpPath, toCpCol_);	
	cpFloatMulArrayMap.insert(std::pair<const std::string, CpMulArray<float> * > (key, arraydata));
	
	arraydata->print();	
	return 0;
}


// ========== UPDATE, WRITE, READ CALLS ========== // 

int CP::updateCp(){
	assert ( cpCommitted == true );
	if(intPod.size()!=0){
		std::map<std::string, const int * const>::iterator it = intPod.begin();
		std::map<std::string, int>::iterator itAsync = intPodAsync.begin();
		for(it = intPod.begin() ; it != intPod.end() ; ++it, ++itAsync){
			itAsync->second = *it->second;
			std::cout << "int: new Async val:" << itAsync->first << ", "<< itAsync->second << std::endl;
		}
	}

	if(doublePod.size()!=0){
		std::map<std::string, const double * const>::iterator it = doublePod.begin();
		std::map<std::string, double>::iterator itAsync = doublePodAsync.begin();
		for(it = doublePod.begin() ; it != doublePod.end() ; ++it, ++itAsync){
			itAsync->second = *it->second;
			std::cout << "dou: new Async val:" << itAsync->first << ", "<< itAsync->second << std::endl;
		}
	}

	if(floatPod.size()!=0){
		std::map<std::string, const float * const >::iterator it = floatPod.begin();
		std::map<std::string, float>::iterator itAsync = floatPodAsync.begin();
		for(it = floatPod.begin() ; it != floatPod.end() ; ++it, ++itAsync){
			itAsync->second = *it->second;
			std::cout << "flo: new Async val:" << itAsync->first << ", "<< itAsync->second << std::endl;
		}
	}

	// ===== GHOST DENSE MAT CALLS ===== // 
	if(cpGhostDenseMat.size()!=0){
		std::map<std::string, CpGhostDenseMat * >::iterator it = cpGhostDenseMat.begin();
		for(it = cpGhostDenseMat.begin(); it != cpGhostDenseMat.end(); ++it){
			it->second->update();
		}
	}
	if(cpGhostDenseMatArray.size()!=0){	
		std::map<const std::string, CpGhostDenseMatArray * >::iterator it = cpGhostDenseMatArray.begin();
		for(it = cpGhostDenseMatArray.begin(); it != cpGhostDenseMatArray.end(); ++it){
			it->second->update();
		}
	}

	// ===== ARRAY CALLS ===== // 
	if(cpIntArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = cpIntArrayMap.begin();
		for (it = cpIntArrayMap.begin(); it != cpIntArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(cpDoubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = cpDoubleArrayMap.begin();
		for (it = cpDoubleArrayMap.begin(); it != cpDoubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(cpFloatArrayMap.size()!=0){
		std::map<std::string, CpArray<float> * >::iterator it = cpFloatArrayMap.begin();
		for (it = cpFloatArrayMap.begin(); it != cpFloatArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(cpIntMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<int> * >::iterator it = cpIntMulArrayMap.begin();
		for (it = cpIntMulArrayMap.begin(); it != cpIntMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(cpDoubleMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<double> * >::iterator it = cpDoubleMulArrayMap.begin();
		for (it = cpDoubleMulArrayMap.begin(); it != cpDoubleMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}
	if(cpFloatMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<float> * >::iterator it = cpFloatMulArrayMap.begin();
		for (it = cpFloatMulArrayMap.begin(); it != cpFloatMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	return 0;
}

int CP::writeCp(){
	assert ( cpCommitted == true );
	int myrank = -1;
	if(cpmpicomm == MPI_COMM_WORLD){
		printf("ERROR: communicator in CP-library is MPI_COMM_WORLD. This do not conform with C/R with automatic FT.\n");	
	}
	MPI_Comm_rank(cpmpicomm, &myrank);
    if(cpPath == "" ){
    		fprintf(stderr, "Error: ckeckpoint path is NULL\n");
		return -1;
    }	


	// ===== GHOST DENSE MAT CALLS ===== // 
	if(cpGhostDenseMat.size()!=0){	
		std::map<std::string, CpGhostDenseMat * >::iterator it = cpGhostDenseMat.begin();
		for(it = cpGhostDenseMat.begin(); it != cpGhostDenseMat.end(); ++it){
			it->second->write();
		}
	}
	if(cpGhostDenseMatArray.size()!=0){	
		std::map<std::string, CpGhostDenseMatArray * >::iterator it = cpGhostDenseMatArray.begin();
		for(it = cpGhostDenseMatArray.begin(); it != cpGhostDenseMatArray.end(); ++it){
			it->second->write();
		}
	}

	// ===== POD CALLS ===== // 
	if(intPod.size() != 0 || doublePod.size() !=0 || floatPod.size() != 0 ){
		char * filename = new char[256];
		sprintf(filename, "%s/POD-rank%d.cp", cpPath.c_str(), myrank);
		FILE *fp1;
		if( NULL == (fp1 = fopen(filename, "w+")) ) {
			fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
			return 0;
		}
		if(intPod.size() != 0){
			std::map<std::string, int>::iterator itAsync = intPodAsync.begin();
			for(itAsync = intPodAsync.begin(); itAsync != intPodAsync.end(); ++itAsync){
				fprintf(fp1, "%s %d \n", itAsync->first.c_str(), itAsync->second);	
			}
		}
		if(doublePod.size() != 0){
			std::map<std::string, double>::iterator itAsync = doublePodAsync.begin();
			for(itAsync = doublePodAsync.begin(); itAsync != doublePodAsync.end(); ++itAsync){
				fprintf(fp1, "%s %.100f \n", itAsync->first.c_str(), itAsync->second);	
			}
		}
		if(floatPod.size() != 0){
			std::map<std::string, float>::iterator itAsync = floatPodAsync.begin();
			for(itAsync = floatPodAsync.begin(); itAsync != floatPodAsync.end(); ++itAsync){
				fprintf(fp1, "%s %.100f \n", itAsync->first.c_str(), itAsync->second);	
			}
		}
		fclose(fp1);
	}


	// ===== ARRAY CALLS ===== // 
	if(cpIntArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = cpIntArrayMap.begin();
		for (it = cpIntArrayMap.begin(); it != cpIntArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(cpDoubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = cpDoubleArrayMap.begin();
		for (it = cpDoubleArrayMap.begin(); it != cpDoubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(cpFloatArrayMap.size()!=0){
		std::map<std::string, CpArray<float> * >::iterator it = cpFloatArrayMap.begin();
		for (it = cpFloatArrayMap.begin(); it != cpFloatArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}

	if(cpIntMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<int> * >::iterator it = cpIntMulArrayMap.begin();
		for (it = cpIntMulArrayMap.begin(); it != cpIntMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(cpDoubleMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<double> * >::iterator it = cpDoubleMulArrayMap.begin();
		for (it = cpDoubleMulArrayMap.begin(); it != cpDoubleMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}	
	if(cpFloatMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<float> * >::iterator it = cpFloatMulArrayMap.begin();
		for (it = cpFloatMulArrayMap.begin(); it != cpFloatMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}

	return 0;
}

int CP::readCp(){
	assert( cpCommitted == true );
	int myrank = 0;
	MPI_Comm_rank(cpmpicomm, &myrank); 
	
	// ===== GHOST DENSE MAT CALLS ===== // 
	if(cpGhostDenseMat.size()!=0){	
		std::map<std::string, CpGhostDenseMat * >::iterator it = cpGhostDenseMat.begin();
		for(it = cpGhostDenseMat.begin(); it != cpGhostDenseMat.end(); ++it){
			cout << "reading densemat " << it->first << endl;
			it->second->read();
		}
	}

	if(cpGhostDenseMatArray.size()!=0){	
		std::map<std::string, CpGhostDenseMatArray * >::iterator it = cpGhostDenseMatArray.begin();
		for(it = cpGhostDenseMatArray.begin(); it != cpGhostDenseMatArray.end(); ++it){
			cout << "reading densematArray " << it->first << endl;
			it->second->read();
		}
	}

	// ===== POD CALLS ===== // 
	if(intPod.size() != 0 || doublePod.size() !=0 || floatPod.size() != 0 ){
			char * cpFile = new char[256];
			sprintf(cpFile, "%s/POD-rank%d.cp", cpPath.c_str(), myrank);
			FILE *fp1;
			if( NULL == (fp1 = fopen(cpFile, "r")) ) {
				fprintf(stderr, "Error: Unable to open file (%s)\n", cpFile);
				return 0;
			}
			if(intPod.size() != 0){
				std::map<const std::string, int>::iterator itAsync = intPodAsync.begin();
				std::map<const std::string, const int * const>::iterator it = intPod.begin();
				for(itAsync = intPodAsync.begin(); itAsync != intPodAsync.end(); ++itAsync, ++it){
					char * tmp1 = new char[256];
					char * tmp2 = new char[256];
					fscanf(fp1, "%s %s", tmp1 , tmp2);
					itAsync->second = atoi(tmp2);	
//					*it->second = itAsync->second;	// TODO: check if it is needed to be read or not TODO: a separeate function could do the trick.
					printf("%d: read int data is: %s %d \n", myrank, itAsync->first.c_str(), itAsync->second);
//					printf("it read data is: %s %d \n", it->first.c_str(), *it->second);
				}
			}
			
			if(doublePod.size() != 0){
				std::map<const std::string, double>::iterator itAsync = doublePodAsync.begin();
				std::map<const std::string, const double * const >::iterator it = doublePod.begin();
				for(itAsync = doublePodAsync.begin(); itAsync != doublePodAsync.end(); ++itAsync, ++it){
					char * tmp1 = new char[256];
					char * tmp2 = new char[256];
					fscanf(fp1, "%s %s", tmp1, tmp2);
					itAsync->second = atof(tmp1);
//					*it->second = itAsync->second;	// TODO: check if it is needed to be read or not	
					printf("%d: read data is: %s %f\n", myrank, itAsync->first.c_str(), itAsync->second);
//					printf("%d: read d data is it: %s %f \n", myrank, it->first.c_str(), *it->second);
				}
			}
			if(floatPod.size() != 0){
				std::map<const std::string, float>::iterator itAsync = floatPodAsync.begin();
				std::map<const std::string, const float * const >::iterator it = floatPod.begin();
				for(itAsync = floatPodAsync.begin(); itAsync != floatPodAsync.end(); ++itAsync, ++it){
					char * tmp1 = new char[256];
					char * tmp2 = new char[256];
				   	fscanf(fp1, "%s %s", tmp1, tmp2);
					itAsync->second = atof(tmp2);
//					*it->second = itAsync->second;	// TODO: check if it is needed to be read or not
					printf("read f data is: %s %f \n", itAsync->first.c_str(), itAsync->second);
				}
			}
			fclose(fp1);
	}
	
	if(cpIntArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = cpIntArrayMap.begin();
		for (it = cpIntArrayMap.begin(); it != cpIntArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
		}
	}
	if(cpDoubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = cpDoubleArrayMap.begin();
		for (it = cpDoubleArrayMap.begin(); it != cpDoubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
		}
	}
	if(cpFloatArrayMap.size()!=0){
		std::map<std::string, CpArray<float> * >::iterator it = cpFloatArrayMap.begin();
		for (it = cpFloatArrayMap.begin(); it != cpFloatArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
		}
	}

	if(cpIntMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<int> * >::iterator it = cpIntMulArrayMap.begin();
		for (it = cpIntMulArrayMap.begin(); it != cpIntMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
			it->second->print();
		}
	}
	if(cpDoubleMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<double> * >::iterator it = cpDoubleMulArrayMap.begin();
		for (it = cpDoubleMulArrayMap.begin(); it != cpDoubleMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
			it->second->print();
		}
	}	
	if(cpFloatMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<float> * >::iterator it = cpFloatMulArrayMap.begin();
		for (it = cpFloatMulArrayMap.begin(); it != cpFloatMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
			it->second->print();
		}
	}

	return 0;
}


	/*if(doubleArray.size()!=0){
		std::map<std::string, std::map<double *, size_t> >::iterator it = doubleArray.begin();
	   	std::map<std::string, std::map<double *, size_t> >::iterator itAsync = doubleArrayAsync.begin();

		std::map<double *, size_t>::iterator itData = doubleArrayData.begin();
		std::map<double *, size_t>::iterator itDataAsync = doubleArrayDataAsync.begin();	  

		for(it = doubleArray.begin(); it != doubleArray.end(); ++it, ++itAsync, ++itData, ++itDataAsync){
	
   		   size_t numElems = itData->second;

//			for(int i = 0; i<numElems; ++i){ 
//				printf("it->first:%s, itData->first[%d]: %f\n" , it->first.c_str(), i, itData->first[i]);
//			}
//
			std::copy(&itData->first[0], &itData->first[numElems], itDataAsync->first);
		
//			for(int i = 0; i<numElems; ++i){ 
//				printf("itDataAsync->first[%d]: %f\n" , i, itDataAsync->first[i]);
//			}

		}

	}*/

/*if(doubleArray.size()!=0){
		std::map<std::string, std::map<double *, size_t> >::iterator itAsync = doubleArrayAsync.begin();
	   	std::map<double *, size_t>::iterator itDataAsync = doubleArrayDataAsync.begin();
		for(itAsync = doubleArrayAsync.begin(); itAsync != doubleArrayAsync.end(); ++itAsync, ++itDataAsync){
	     	char *filename = new char[256];	
			sprintf(filename, "%s/%s-rank%d.cp", cpPath, itAsync->first.c_str(), myrank);	
			FILE *fp;
			if( NULL == (fp = fopen(filename, "w+")) ) {
				fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
				return -1;
			}
			for(int i = 0 ; i < itDataAsync->second; ++i) {
//				fprintf(fp, "%f\n", itDataAsync->first[i]);	
//				printf("%f\n", itDataAsync->first[i]);	
			}	

			fwrite( &itDataAsync->second, sizeof(size_t), 1, fp );
			fwrite( itDataAsync->first, sizeof(double), itDataAsync->second, fp );
			fclose(fp);	
		}
	}*/
/*

	if(doubleArray.size()!=0){

		std::map<std::string, std::map<double *, size_t> >::iterator itAsync = doubleArrayAsync.begin();
	   	std::map<double *, size_t>::iterator itData = doubleArrayData.begin();
	   	std::map<double *, size_t>::iterator itDataAsync = doubleArrayDataAsync.begin();
		
		for(itAsync = doubleArrayAsync.begin(); itAsync != doubleArrayAsync.end(); ++itAsync, ++itData, ++itDataAsync){
	     	char *filename = new char[256];	
			sprintf(filename, "%s/%s-rank%d.cp", cpPath, itAsync->first.c_str(), myrank);	
			FILE *fp;
			if( NULL == (fp = fopen(filename, "r")) ) {
				fprintf(stderr, "Error: Unable to open file (%s)\n", filename);
				return -1;
			}
			
			size_t tempArrayLength = 0; 
			fread( &tempArrayLength, sizeof(size_t), 1, fp );
			fread( &itDataAsync->first[0], sizeof(double), itDataAsync->second, fp );
			fclose(fp);
			assert (tempArrayLength == itDataAsync->second);

			std::copy(&itDataAsync->first[0], &itDataAsync->first[tempArrayLength], itData->first);	

//			for(int i = 0 ; i < itDataAsync->second; ++i) {
//				printf("%d: %s after reading: %f \n", myrank, itAsync->first.c_str(), itData->first[i]);	
//			}	
		
		}
	}	     */




