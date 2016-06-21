#ifndef __CPLIB_H__
#define __CPLIB_H__

#include <map>
#include <vector>
#include <cstdlib>
#include <string>
#include <ghost.h>
#include <ghost/types.h>
#include <iostream>
#include <complex>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "cp_array.h"
#include "cp_ghostdensemat.h"
#include "cpPOD.h"
#include <mpi.h>

#ifdef SCR
extern "C"{
	#include <scr.h>
}
#endif

using namespace std;
using std::complex;
class CP
{
protected:	
	MPI_Comm cpmpicomm;
	std::string cpPath;
	bool useSCR;							// true = enabled, false = disabled 
	bool cpCommitted;

private:
	// ===== POD ===== 
	std::map<const std::string, CpPod<int> * > intPodMap;
	std::map<const std::string, CpPod<double> *> doublePodMap;
	std::map<const std::string, CpPod<float> *> floatPodMap;
	std::map<const std::string, CpPod< complex<double> >* > doubleComplexPodMap;
	
	// ===== GHOST-VECTORS ===== 
	std::map<const std::string, CpGhostDenseMat *> ghostDenseMatMap;
	std::map<const std::string, CpGhostDenseMatArray *> ghostDenseMatArrayMap;

	// ===== ARRAYS ===== 
 	std::map<const std::string, CpArray<int> * > 	intArrayMap;
	std::map<const std::string, CpArray<double> * > doubleArrayMap;
	std::map<const std::string, CpArray<float> * > 	floatArrayMap;
	std::map<const std::string, CpArray<complex<double> > * > doubleComplexArrayMap;

	// ===== MULTI ARRAYS ===== 
	std::map<const std::string, CpMulArray<int> * > 	intMulArrayMap;
	std::map<const std::string, CpMulArray<double> * > 	doubleMulArrayMap;
	std::map<const std::string, CpMulArray<float> * > 	floatMulArrayMap;


	int CP_ADD_POD_INT		(const std::string key, int * const value);
	int CP_ADD_POD_DOUBLE	(const std::string key, double * const value); 
	int CP_ADD_POD_FLOAT	(const std::string key, float * const value); 
	
	int CP_ADD_POD_DOUBLE_COMPLEX	(const std::string key, complex<double> * const value); 

	int CP_ADD_ARRAY_INT	(const std::string key, int * const val_array, const size_t array_size );
	int CP_ADD_ARRAY_DOUBLE	(const std::string key, double * const val_array, const size_t array_size );
	int CP_ADD_ARRAY_FLOAT	(const std::string key, float * const val_array, const size_t array_size );

	int CP_ADD_ARRAY_DOUBLE_COMPLEX	(const std::string key, complex<double>  * const val_array, const size_t array_size );


	int CP_ADD_MULTIARRAY_INT	(const std::string key, int** const ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );
	int CP_ADD_MULTIARRAY_DOUBLE(const std::string key, double** const ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );
	int CP_ADD_MULTIARRAY_FLOAT(const std::string key, float** const ptr, const size_t n_rows, const size_t n_cols, const int toCpCol_ );

public:
	CP();
	~CP();
	int setCpPath(const std::string cpPath_);
	int setComm(const MPI_Comm FT_Comm);
	void enableSCR();
	void commit();

	int CP_ADD_GHOST_DENSEMAT(const std::string key, ghost_densemat * const value);
	int CP_ADD_GHOST_DENSEMAT_ARRAY(const std::string key, ghost_densemat ** const value, int num_vecs, const int toCpVec_ = ALL);
	
	int CP_ADD_POD(const std::string key, int * const val_ptr);
	int CP_ADD_POD(const std::string key, double * const val_ptr);
	int CP_ADD_POD(const std::string key, complex<double>  * const val_ptr);
	int CP_ADD_POD(const std::string key, float * const val_ptr);

	int CP_ADD_ARRAY(const std::string key,  int 	* const val_array, const size_t array_size);
	int CP_ADD_ARRAY(const std::string key,  double * const val_array, const size_t array_size);
	int CP_ADD_ARRAY(const std::string key,  float * const val_array, const size_t array_size);
	
	int CP_ADD_ARRAY(const std::string key,  complex<double>  * const val_array, const size_t array_size);

	int CP_ADD_MULTIARRAY(const std::string key, int** const ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);
	int CP_ADD_MULTIARRAY(const std::string key, double** const ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);
	int CP_ADD_MULTIARRAY(const std::string key, float** const ptr,  const size_t n_rows,  const size_t n_cols, const int toCpCol_ = ALL);

	int update();	// this function should update the values of CP
	int write();		// this function should write the updated CP
	int read();

	int print();
};

#endif

