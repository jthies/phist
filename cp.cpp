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

int CP::CP_ADD_VEC(std::string key, ghost_densemat * value){
	assert (cpCommitted == false)	;
	SAFE_INSERT( vec.insert( std::pair<std::string, ghost_densemat *>(key, value)) );

	// === make async copy of vector === 
	ghost_densemat_create(&vector_async, value->context, value->traits);
	ghost_densemat_init_densemat(vector_async, value, 0, 0);
	SAFE_INSERT( vecAsync.insert( std::pair<std::string, ghost_densemat * > (key, vector_async)) );
	return 0;
}


int CP::CP_ADD_MULTI_VEC(const std::string key, ghost_densemat ** value, int num_vecs){
	assert (cpCommitted == false);

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

int CP::CP_ADD_ARRAY_INT(const std::string key, const int * const val_array, const size_t array_size ){
	printf("Type of array is int\n");

	CpArray<int> * arraydata = new CpArray<int>[1];
	arraydata->add(key, val_array, array_size, cpmpicomm, cpPath );	
	arraydata->print();	
	CpIntArrayMap.insert(std::pair<const std::string, CpArray<int> * >(key, arraydata));	

/*     SAFE_INSERT( doubleArrayData.insert( std::pair<double *, size_t> (val_array, array_size)) );
     SAFE_INSERT( doubleArray.insert( std::pair<std::string, std::map<double *, size_t> > (key, doubleArrayData)) );

	size_t num_elems = array_size;
	double * valArrayAsync = new double [num_elems];
	std::copy(&val_array[0], &val_array[num_elems], valArrayAsync);

	SAFE_INSERT( doubleArrayDataAsync.insert( std::pair<double *, size_t> (valArrayAsync, array_size)) );
     SAFE_INSERT( doubleArrayAsync.insert( std::pair<std::string, std::map<double * , size_t > > (key, doubleArrayDataAsync)) );
*/
	return 0;	
}

int CP::CP_ADD_ARRAY_DOUBLE(const std::string key, const double * const val_array, const size_t array_size ){
	printf("Type of array is double\n");
	CpArray<double> * arraydata = new CpArray<double>[1];
	arraydata->add(key, val_array, array_size , cpmpicomm, cpPath);
//	arraydata->print();
	CpDoubleArrayMap.insert(std::pair<const std::string, CpArray<double> * >(key, arraydata));	


    /* SAFE_INSERT( doubleArrayData.insert( std::pair<double *, size_t> (val_array, array_size)) );
     SAFE_INSERT( doubleArray.insert( std::pair<std::string, std::map<double *, size_t> > (key, doubleArrayData)) );

	size_t num_elems = array_size;
	double * valArrayAsync = new double [num_elems];
	std::copy(&val_array[0], &val_array[num_elems], valArrayAsync);

	SAFE_INSERT( doubleArrayDataAsync.insert( std::pair<double *, size_t> (valArrayAsync, array_size)) );
     SAFE_INSERT( doubleArrayAsync.insert( std::pair<std::string, std::map<double * , size_t > > (key, doubleArrayDataAsync)) );
*/
	return 0;	
}



// ========== MULTI ARRAY CALLS ========== //

int CP::CP_ADD_MULTIARRAY(const std::string key, int ** ptr, size_t nRows, size_t n_cols){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_INT( key,  ptr, nRows, n_cols);
	return 0;
}

int CP::CP_ADD_MULTIARRAY_INT(const std::string key, int ** ptr, size_t nRows, size_t n_cols){
	assert (cpCommitted == false);
	printf ("Type of array is int\n");	
	
	//CpMulArray<int> * arraydata = new CpArray<int>[1];
//	arraydata->ADD(ptr, nRows, n_cols);	
//	arraydata->Print(ptr);	
	return 0;
}

int CP::CP_ADD_MULTIARRAY(const std::string key, double ** ptr, size_t nRows, size_t n_cols){
	assert (cpCommitted == false);
	CP_ADD_MULTIARRAY_DOUBLE( key,  ptr, nRows, n_cols);
	return 0;
}

int CP::CP_ADD_MULTIARRAY_DOUBLE(const std::string key, double ** ptr, size_t nRows, size_t n_cols){
	assert (cpCommitted == false);
	printf ("Type of array is double\n");	
	
//	CpArray<double> * arraydata = new CpArray<double>[1];
//	arraydata->ADD(ptr, nRows, n_cols);	
//	arraydata->Print(ptr);	
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

	if(vec.size()!=0){
		std::map<std::string, ghost_densemat *>::iterator it = vec.begin();
		std::map<std::string, ghost_densemat *>::iterator itAsync = vecAsync.begin();
		for(it = vec.begin() ; it != vec.end() ; ++it, ++itAsync){
	  		ghost_densemat_init_densemat(itAsync->second, it->second, 0,0); 
			usleep(100000);
			std::cout << "vec: vec is updated" << itAsync->first << std::endl;
//			char * vecstr;
//			ghost_densemat_string(&vecstr, itAsync->second);
//			printf("y: \n%s\n", vecstr);	
		}
	}

	if(CpIntArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = CpIntArrayMap.begin();
		for (it = CpIntArrayMap.begin(); it != CpIntArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->update();
		}
	}

	if(CpDoubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = CpDoubleArrayMap.begin();
		for (it = CpDoubleArrayMap.begin(); it != CpDoubleArrayMap.end(); ++it){
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
	if(vecAsync.size()!=0){	
		std::map<std::string, ghost_densemat *>::iterator itVecAsync = vecAsync.begin();
		for(itVecAsync = vecAsync.begin(); itVecAsync != vecAsync.end(); ++itVecAsync){
			char * filename = new char[256];
			sprintf(filename, "%s/%s-rank%d.cp", cpPath, itVecAsync->first.c_str(), myrank);
			printf("filename: %s\n", filename);
			itVecAsync->second->toFile( itVecAsync->second, filename, 0);
		}
	}
	if(intPod.size() != 0 || doublePod.size() !=0 || floatPod.size() != 0 ){
		char * filename = new char[256];
		sprintf(filename, "%s/POD-rank%d.cp", cpPath, myrank);
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
	
	if(CpIntArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = CpIntArrayMap.begin();
		for (it = CpIntArrayMap.begin(); it != CpIntArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->write();
		}
	}
	if(CpDoubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = CpDoubleArrayMap.begin();
		for (it = CpDoubleArrayMap.begin(); it != CpDoubleArrayMap.end(); ++it){
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
	if(vecAsync.size()!=0){	
		std::map<std::string, ghost_densemat *>::iterator itVecAsync = vecAsync.begin();
		std::map<std::string, ghost_densemat *>::iterator itVec = vec.begin();
		for(itVecAsync = vecAsync.begin(); itVecAsync != vecAsync.end(); ++itVecAsync, ++itVec){
			char * filename = new char[256];
			sprintf(filename, "%s/%s-rank%d.cp", cpPath, itVecAsync->first.c_str(), myrank);
			printf("%d: filename: %s\n", myrank, filename);
			itVecAsync->second->fromFile( itVecAsync->second, filename, 0);
	  		ghost_densemat_init_densemat(itVec->second, itVecAsync->second, 0,0);
			char * readstr = new char[256];	
			ghost_densemat_string (&readstr, itVec->second);
//			printf("itVec->second: \n%s\n", readstr);
		}
	}
	if(intPod.size() != 0 || doublePod.size() !=0 || floatPod.size() != 0 ){
			char * cpFile = new char[256];
			sprintf(cpFile, "%s/POD-rank%d.cp", cpPath, myrank);
			FILE *fp1;
			if( NULL == (fp1 = fopen(cpFile, "r")) ) {
				fprintf(stderr, "Error: Unable to open file (%s)\n", cpFile);
				return 0;
			}
			if(intPod.size() != 0){
				std::map<const std::string, int>::iterator itAsync = intPodAsync.begin();
				std::map<const std::string, const int * const>::iterator it = intPod.begin();
				for(itAsync = intPodAsync.begin(); itAsync != intPodAsync.end(); ++itAsync, ++it){
					char * tmp = new char[256];
					fscanf(fp1, "%s %s", itAsync->first, tmp);
					itAsync->second = atoi(tmp);	
//					*it->second = itAsync->second;	// TODO: check if it is needed to be read or not
					printf("%d: read int data is: %s %d \n", myrank, itAsync->first.c_str(), itAsync->second);
//					printf("it read data is: %s %d \n", it->first.c_str(), *it->second);
				}
			}
			
			if(doublePod.size() != 0){
				std::map<const std::string, double>::iterator itAsync = doublePodAsync.begin();
				std::map<const std::string, const double * const >::iterator it = doublePod.begin();
				for(itAsync = doublePodAsync.begin(); itAsync != doublePodAsync.end(); ++itAsync, ++it){
					char * tmp = new char[256];
					fscanf(fp1, "%s %s", itAsync->first, tmp);
					itAsync->second = atof(tmp);
//					*it->second = itAsync->second;	// TODO: check if it is needed to be read or not	
					printf("%d: read data is: %s %f\n", myrank, itAsync->first.c_str(), itAsync->second);
//					printf("%d: read d data is it: %s %f \n", myrank, it->first.c_str(), *it->second);
				}
			}
			if(floatPod.size() != 0){
				std::map<const std::string, float>::iterator itAsync = floatPodAsync.begin();
				std::map<const std::string, const float * const >::iterator it = floatPod.begin();
				for(itAsync = floatPodAsync.begin(); itAsync != floatPodAsync.end(); ++itAsync, ++it){
					char * tmp = new char[256];
				   	fscanf(fp1, "%s %f", itAsync->first, tmp);
					itAsync->second = atof(tmp);
//					*it->second = itAsync->second;	// TODO: check if it is needed to be read or not
					printf("read f data is: %s %f \n", itAsync->first.c_str(), itAsync->second);
				}
			}
			fclose(fp1);
	}
	
	if(CpIntArrayMap.size()!=0){
		std::map<std::string, CpArray<int> * >::iterator it = CpIntArrayMap.begin();
		for (it = CpIntArrayMap.begin(); it != CpIntArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
		}
	}
	if(CpDoubleArrayMap.size()!=0){
		std::map<std::string, CpArray<double> * >::iterator it = CpDoubleArrayMap.begin();
		for (it = CpDoubleArrayMap.begin(); it != CpDoubleArrayMap.end(); ++it){
			cout << it->first << endl;
			it->second->read();
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




