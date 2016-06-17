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
#ifdef SCR
	printf(	"scr is defined in CP.cpp\n");
#endif
#ifndef SCR
	printf(	"scr is NOT defined in CP.cpp\n");
#endif


}

CP::~CP(){
	printf("============= SCR_Finalize is done1 ============= \n");
#ifdef SCR
	printf("============= SCR_Finalize is done ============= \n");
	SCR_Finalize();
#endif
}

int CP::setCpPath(const std::string cpPath_){
	cpPath = cpPath_;
	return 0;
}
	
int CP::setComm(const MPI_Comm FT_Comm){
	cpmpicomm = FT_Comm;
#ifdef SCR
	SCR_Init(&cpmpicomm);
	printf("============= SCR_init is done with cpmpicomm ============= \n");
#endif
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
	ghostDenseMatMap.insert(std::pair<const std::string, CpGhostDenseMat * >(key, cpdata));	
	return 0;
}


int CP::CP_ADD_GHOST_DENSEMAT_ARRAY(const std::string key, ghost_densemat ** const value, const int nVecs, const int toCpVec_)
{
	assert (cpCommitted == false);
	CpGhostDenseMatArray * cpdata= new CpGhostDenseMatArray[1];
	cpdata->add(key, value , nVecs, cpmpicomm, cpPath, toCpVec_);
	ghostDenseMatArrayMap.insert( std::pair< const std::string, CpGhostDenseMatArray * >(key, cpdata));
	return 0;
}
	

// ========== POD CALLS ========== //

int CP::CP_ADD_POD(const std::string key, int * const val_ptr){
	assert (cpCommitted == false);
	CP_ADD_POD_INT(key, val_ptr);
	return 0;
}

int CP::CP_ADD_POD(const std::string key, double * const val_ptr){
	assert (cpCommitted == false);
	CP_ADD_POD_DOUBLE(key, val_ptr);	
	return 0;
}

int CP::CP_ADD_POD(const std::string key, complex<double> * const val_ptr){
	assert (cpCommitted == false);
	CP_ADD_POD_DOUBLE_COMPLEX(key, val_ptr);	
	return 0;
}

int CP::CP_ADD_POD(const std::string key, float * const val_ptr){
	assert (cpCommitted == false)	;
	CP_ADD_POD_FLOAT(key, val_ptr);	
	return 0;
}

int CP::CP_ADD_POD_INT(const std::string key, int * const value){

	CpPod<int> * podData = new CpPod<int>[1];
	podData->add(key, value, cpmpicomm, cpPath );		
	podData->print();	
	intPodMap.insert(std::pair<const std::string, CpPod<int> * >(key, podData));	
	
	return 0;
}

int CP::CP_ADD_POD_DOUBLE(const std::string key, double * const value){	 
	CpPod<double> * podData = new CpPod<double>[1];
	podData->add(key, value, cpmpicomm, cpPath );		
	podData->print();	
	doublePodMap.insert(std::pair<const std::string, CpPod<double> * >(key, podData));	
	return 0;
}

int CP::CP_ADD_POD_DOUBLE_COMPLEX(const std::string key, complex<double> * const value){	 
	CpPod<complex<double> > * podData = new CpPod<complex<double> >[1];
	podData->add(key, value, cpmpicomm, cpPath );		
	podData->print();	
	doubleComplexPodMap.insert(std::pair<const std::string, CpPod<complex<double> > * >(key, podData));	
	return 0;
}

int CP::CP_ADD_POD_FLOAT(const std::string key, float * const value){ 
	CpPod<float> * podData = new CpPod<float>[1];
	podData->add(key, value, cpmpicomm, cpPath );		
	podData->print();	
	floatPodMap.insert(std::pair<const std::string, CpPod<float> * >(key, podData));	
	return 0;
}


// ========== SINGLE ARRAY CALLS ========== //

int CP::CP_ADD_ARRAY(const std::string key, int * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_INT(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY(const std::string key, double * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_DOUBLE(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY(const std::string key, float * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_FLOAT(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY(const std::string key, complex<double>  * const val_array, const size_t array_size){
	assert (cpCommitted == false);
	CP_ADD_ARRAY_DOUBLE_COMPLEX(key, val_array, array_size);
	return 0;
}

int CP::CP_ADD_ARRAY_INT(const std::string key, int * const val_array, const size_t array_size ){
	printf("Type of array is int\n");

	CpArray<int> * arraydata = new CpArray<int>[1];
	arraydata->add(key, val_array, array_size, cpmpicomm, cpPath );		
	arraydata->print();	
	intArrayMap.insert(std::pair<const std::string, CpArray<int> * >(key, arraydata));	

	return 0;	
}

int CP::CP_ADD_ARRAY_DOUBLE(const std::string key, double * const val_array, const size_t array_size ){
	printf("Type of array is double\n");
	CpArray<double> * arraydata = new CpArray<double>[1];
	arraydata->add(key, val_array, array_size , cpmpicomm, cpPath);
//	arraydata->print();
	doubleArrayMap.insert(std::pair<const std::string, CpArray<double> * >(key, arraydata));	

	return 0;	
}

int CP::CP_ADD_ARRAY_FLOAT(const std::string key, float * const val_array, const size_t array_size ){
	printf("Type of array is float\n");
	CpArray<float> * arraydata = new CpArray<float>[1];
	arraydata->add(key, val_array, array_size , cpmpicomm, cpPath);
//	arraydata->print();
	floatArrayMap.insert(std::pair<const std::string, CpArray<float> * >(key, arraydata));	

	return 0;	
}

int CP::CP_ADD_ARRAY_DOUBLE_COMPLEX	(const std::string key, complex<double>  * const val_array, const size_t array_size ){
	printf("Type of array is complex<double> ");
	CpArray<complex<double>  > * arraydata = new CpArray<complex<double>  >[1];
	arraydata->add(key, val_array, array_size , cpmpicomm, cpPath);
//	arraydata->print();
	doubleComplexArrayMap.insert(std::pair<const std::string, CpArray<complex<double>  > * >(key, arraydata));	

	return 0;
}



// ========== MULTI ARRAY CALLS ========== //

int CP::CP_ADD_MULTIARRAY(const std::string key, int ** const ptr, const size_t nRows, const size_t nCols, const int toCpCol_ ){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_INT( key,  ptr, nRows, nCols, toCpCol_);
	return 0;
}

int CP::CP_ADD_MULTIARRAY(const std::string key, double** const ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_DOUBLE( key,  ptr, nRows, nCols, toCpCol_);
	return 0;
}

int CP::CP_ADD_MULTIARRAY(const std::string key, float** const ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_FLOAT( key,  ptr, nRows, nCols, toCpCol_);
	return 0;
}

int CP::CP_ADD_MULTIARRAY_INT(const std::string key, int** const ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	
	CpMulArray<int> * arraydata = new CpMulArray<int>[1];
	arraydata->add(key, ptr, nRows, nCols, cpmpicomm, cpPath, toCpCol_);	
	intMulArrayMap.insert(std::pair<const std::string, CpMulArray<int> * > (key, arraydata));

	arraydata->print();
	return 0;
}

int CP::CP_ADD_MULTIARRAY_DOUBLE(const std::string key, double ** const ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	printf ("Adding DOUBLE MULTIARRAY\n");	

	CpMulArray<double> * arraydata = new CpMulArray<double>[1];
	arraydata->add(key, ptr, nRows, nCols, cpmpicomm, cpPath, toCpCol_);	
	doubleMulArrayMap.insert(std::pair<const std::string, CpMulArray<double> * > (key, arraydata));
	
	arraydata->print();	
	return 0;
}


int CP::CP_ADD_MULTIARRAY_FLOAT(const std::string key, float ** const ptr, const size_t nRows, const size_t nCols, const int toCpCol_){
	printf ("Adding FLOAT MULTIARRAY\n");	

	CpMulArray<float> * arraydata = new CpMulArray<float>[1];
	arraydata->add(key, ptr, nRows, nCols, cpmpicomm, cpPath, toCpCol_);	
	floatMulArrayMap.insert(std::pair<const std::string, CpMulArray<float> * > (key, arraydata));
	
	arraydata->print();	
	return 0;
}


// ========== UPDATE, WRITE, READ CALLS ========== // 

int CP::update(){
	assert ( cpCommitted == true );
	// ===== POD CALLS ===== //
	if(intPodMap.size()!=0){
		std::map<const std::string, CpPod<int> * >::iterator it = intPodMap.begin();	
		for(it = intPodMap.begin() ; it != intPodMap.end() ; ++it){
			it->second->update();
		}
	}

	if(doublePodMap.size()!=0){
		std::map<const std::string, CpPod<double> * >::iterator it = doublePodMap.begin();	
		for(it = doublePodMap.begin() ; it != doublePodMap.end() ; ++it){
			it->second->update();
		}
	}

	if(floatPodMap.size()!=0){
		std::map<const std::string, CpPod<float> * >::iterator it = floatPodMap.begin();	
		for(it = floatPodMap.begin() ; it != floatPodMap.end() ; ++it){
			it->second->update();
		}
	}

	if(doubleComplexPodMap.size()!=0){
		std::map<const std::string, CpPod<complex<double>  > * >::iterator it = doubleComplexPodMap.begin();	
		for(it = doubleComplexPodMap.begin() ; it != doubleComplexPodMap.end() ; ++it){
			it->second->update();
		}
	}

	// ===== GHOST DENSE MAT CALLS ===== // 
	if(ghostDenseMatMap.size()!=0){
		std::map<const std::string, CpGhostDenseMat * >::iterator it = ghostDenseMatMap.begin();
		for(it = ghostDenseMatMap.begin(); it != ghostDenseMatMap.end(); ++it){
			it->second->update();
		}
	}
	if(ghostDenseMatArrayMap.size()!=0){	
		std::map<const std::string, CpGhostDenseMatArray * >::iterator it = ghostDenseMatArrayMap.begin();
		for(it = ghostDenseMatArrayMap.begin(); it != ghostDenseMatArrayMap.end(); ++it){
			it->second->update();
		}
	}

	// ===== ARRAY CALLS ===== // 
	if(intArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = intArrayMap.begin();
		for (it = intArrayMap.begin(); it != intArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(doubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = doubleArrayMap.begin();
		for (it = doubleArrayMap.begin(); it != doubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(floatArrayMap.size()!=0){
		std::map<std::string, CpArray<float> * >::iterator it = floatArrayMap.begin();
		for (it = floatArrayMap.begin(); it != floatArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(intMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<int> * >::iterator it = intMulArrayMap.begin();
		for (it = intMulArrayMap.begin(); it != intMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(doubleMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<double> * >::iterator it = doubleMulArrayMap.begin();
		for (it = doubleMulArrayMap.begin(); it != doubleMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}
	if(floatMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<float> * >::iterator it = floatMulArrayMap.begin();
		for (it = floatMulArrayMap.begin(); it != floatMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	return 0;
}

int CP::write(){
#ifdef SCR
		int valid = 0;
		SCR_Start_checkpoint();
		printf("SCR_Start_checkpoint is called \n");
#endif
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
	if(ghostDenseMatMap.size()!=0){	
		std::map<const std::string, CpGhostDenseMat * >::iterator it = ghostDenseMatMap.begin();
		for(it = ghostDenseMatMap.begin(); it != ghostDenseMatMap.end(); ++it){
			it->second->write();
		}
	}
	if(ghostDenseMatArrayMap.size()!=0){	
		std::map<std::string, CpGhostDenseMatArray * >::iterator it = ghostDenseMatArrayMap.begin();
		for(it = ghostDenseMatArrayMap.begin(); it != ghostDenseMatArrayMap.end(); ++it){
			it->second->write();
		}
	}

	// ===== POD CALLS ===== //
	if(intPodMap.size()!=0){
		std::map<const std::string, CpPod<int> * >::iterator it = intPodMap.begin();	
		for(it = intPodMap.begin() ; it != intPodMap.end() ; ++it){
			it->second->write();
		}
	}

	if(doublePodMap.size()!=0){
		std::map<const std::string, CpPod<double> * >::iterator it = doublePodMap.begin();	
		for(it = doublePodMap.begin() ; it != doublePodMap.end() ; ++it){
			it->second->write();
		}
	}

	if(floatPodMap.size()!=0){
		std::map<const std::string, CpPod<float> * >::iterator it = floatPodMap.begin();	
		for(it = floatPodMap.begin() ; it != floatPodMap.end() ; ++it){
			it->second->write();
		}
	}

	if(doubleComplexPodMap.size()!=0){
		std::map<const std::string, CpPod<complex<double>  > * >::iterator it = doubleComplexPodMap.begin();	
		for(it = doubleComplexPodMap.begin() ; it != doubleComplexPodMap.end() ; ++it){
			it->second->write();
		}
	}  

	// ===== ARRAY CALLS ===== // 
	if(intArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = intArrayMap.begin();
		for (it = intArrayMap.begin(); it != intArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(doubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = doubleArrayMap.begin();
		for (it = doubleArrayMap.begin(); it != doubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(floatArrayMap.size()!=0){
		std::map<std::string, CpArray<float> * >::iterator it = floatArrayMap.begin();
		for (it = floatArrayMap.begin(); it != floatArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}

	if(intMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<int> * >::iterator it = intMulArrayMap.begin();
		for (it = intMulArrayMap.begin(); it != intMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(doubleMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<double> * >::iterator it = doubleMulArrayMap.begin();
		for (it = doubleMulArrayMap.begin(); it != doubleMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}	
	if(floatMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<float> * >::iterator it = floatMulArrayMap.begin();
		for (it = floatMulArrayMap.begin(); it != floatMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}

#ifdef SCR
	valid = 1;
	SCR_Complete_checkpoint(valid);			// valid flag should be 1 if every CP is successfully written. TODO: check the flag from each write call.
#endif

	return 0;
}

int CP::read(){
	assert( cpCommitted == true );
	int myrank = 0;
	MPI_Comm_rank(cpmpicomm, &myrank); 
	
	// ===== GHOST DENSE MAT CALLS ===== // 
	if(ghostDenseMatMap.size()!=0){	
		std::map<std::string, CpGhostDenseMat * >::iterator it = ghostDenseMatMap.begin();
		for(it = ghostDenseMatMap.begin(); it != ghostDenseMatMap.end(); ++it){
			it->second->read();
		}
	}

	if(ghostDenseMatArrayMap.size()!=0){	
		std::map<std::string, CpGhostDenseMatArray * >::iterator it = ghostDenseMatArrayMap.begin();
		for(it = ghostDenseMatArrayMap.begin(); it != ghostDenseMatArrayMap.end(); ++it){
			it->second->read();
		}
	}

	// ===== POD CALLS ===== // 
	if(intPodMap.size()!=0){
		std::map<const std::string, CpPod<int> * >::iterator it = intPodMap.begin();	
		for(it = intPodMap.begin() ; it != intPodMap.end() ; ++it){
			it->second->read();
			//it->second->print();
		}
	}

	if(doublePodMap.size()!=0){
		std::map<const std::string, CpPod<double> * >::iterator it = doublePodMap.begin();	
		for(it = doublePodMap.begin() ; it != doublePodMap.end() ; ++it){
			it->second->read();
	//		it->second->print();
		}
	}

	if(floatPodMap.size()!=0){
		std::map<const std::string, CpPod<float> * >::iterator it = floatPodMap.begin();	
		for(it = floatPodMap.begin() ; it != floatPodMap.end() ; ++it){
			it->second->read();
	//		it->second->print();
		}
	}

	if(doubleComplexPodMap.size()!=0){
		std::map<const std::string, CpPod<complex<double>  > * >::iterator it = doubleComplexPodMap.begin();	
		for(it = doubleComplexPodMap.begin() ; it != doubleComplexPodMap.end() ; ++it){
			it->second->read();
	//		it->second->print();
		}
	} 

	// ===== ARRAY CALLS ===== // 
	if(intArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = intArrayMap.begin();
		for (it = intArrayMap.begin(); it != intArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
	//		it->second->print();
		}
	}
	if(doubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = doubleArrayMap.begin();
		for (it = doubleArrayMap.begin(); it != doubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
	//		it->second->print();
		}
	}
	if(floatArrayMap.size()!=0){
		std::map<std::string, CpArray<float> * >::iterator it = floatArrayMap.begin();
		for (it = floatArrayMap.begin(); it != floatArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
	//		it->second->print();
		}
	}

	if(intMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<int> * >::iterator it = intMulArrayMap.begin();
		for (it = intMulArrayMap.begin(); it != intMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
	//		it->second->print();
		}
	}
	if(doubleMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<double> * >::iterator it = doubleMulArrayMap.begin();
		for (it = doubleMulArrayMap.begin(); it != doubleMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
	//		it->second->print();
		}
	}	
	if(floatMulArrayMap.size()!=0){
		std::map<std::string, CpMulArray<float> * >::iterator it = floatMulArrayMap.begin();
		for (it = floatMulArrayMap.begin(); it != floatMulArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
	//		it->second->print();
		}
	}

	return 0;
}


