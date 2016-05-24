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
#include "cp_array.h"
#include "cp_ghostdensemat.h"

using namespace std;

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
	std::map<const std::string, CpGhostDenseMat *> cpGhostDenseMat;
	std::map<const std::string, CpGhostDenseMatArray *> cpGhostDenseMatArray;

	//ghost_densemat *vector_async;

	// ===== ARRAYS ===== 

	std::map<const std::string, CpArray<int> * > 	cpIntArrayMap;
	std::map<const std::string, CpArray<double> * > cpDoubleArrayMap;
	std::map<const std::string, CpArray<float> * > 	cpFloatArrayMap;
	
	std::map<const std::string, CpMulArray<int> * > 	cpIntMulArrayMap;
	std::map<const std::string, CpMulArray<double> * > 	cpDoubleMulArrayMap;
	std::map<const std::string, CpMulArray<float> * > 	cpFloatMulArrayMap;

	int CP_ADD_POD_INT		(const std::string key, const int * const value);
	int CP_ADD_POD_DOUBLE	(const std::string key, const double * const value); 
	int CP_ADD_POD_FLOAT	(const std::string key, const float * const value); 

	int CP_ADD_ARRAY_INT	(const std::string key, const int * const val_array, const size_t array_size );
	int CP_ADD_ARRAY_DOUBLE	(const std::string key, const double * const val_array, const size_t array_size );
	int CP_ADD_ARRAY_FLOAT	(const std::string key, const float * const val_array, const size_t array_size );

	int CP_ADD_MULTIARRAY_INT	(const std::string key, const int* const* ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );
	int CP_ADD_MULTIARRAY_DOUBLE(const std::string key, const double* const* ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );
	int CP_ADD_MULTIARRAY_FLOAT(const std::string key, const float* const* ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );

public:
	CP();
	~CP();
	int setCpPath(const std::string cpPath_);
	int setComm(const MPI_Comm FT_Comm);
	void commit();

	int CP_ADD_GHOST_DENSEMAT(const std::string key, ghost_densemat * const value);
	int CP_ADD_GHOST_DENSEMAT_ARRAY(const std::string key, ghost_densemat ** const value, int num_vecs, const int toCpVec_ = ALL);
	
	int CP_ADD_POD(const std::string key, const int * const val_ptr);
	int CP_ADD_POD(const std::string key, const double * const val_ptr);
	int CP_ADD_POD(const std::string key, const float * const val_ptr);

	int CP_ADD_ARRAY(const std::string key,  const int 	* const val_array, const size_t array_size);
	int CP_ADD_ARRAY(const std::string key,  const double * const val_array, const size_t array_size);
	int CP_ADD_ARRAY(const std::string key,  const float * const val_array, const size_t array_size);

	int CP_ADD_MULTIARRAY(const std::string key, const int* const* ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);
	int CP_ADD_MULTIARRAY(const std::string key, const double* const* ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);
	int CP_ADD_MULTIARRAY(const std::string key, const float* const* ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);

	int updateCp();	// this function should update the values of CP
	int writeCp();		// this function should write the updated CP
	int readCp();

	int print();
};

#endif



