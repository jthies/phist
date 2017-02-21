/*******************************************************************************
 * *  Copyright (C) 2009-2015 Intel Corporation. All Rights Reserved.
 * *  The information and material ("Material") provided below is owned by Intel
 * *  Corporation or its suppliers or licensors, and title to such Material remains
 * *  with Intel Corporation or its suppliers or licensors. The Material contains
 * *  proprietary information of Intel or its suppliers and licensors. The Material
 * *  is protected by worldwide copyright laws and treaty provisions. No part of
 * *  the Material may be copied, reproduced, published, uploaded, posted,
 * *  transmitted, or distributed in any way without Intel's prior express written
 * *  permission. No license under any patent, copyright or other intellectual
 * *  property rights in the Material is granted to or conferred upon you, either
 * *  expressly, by implication, inducement, estoppel or otherwise. Any license
 * *  under such intellectual property rights must be express and approved by Intel
 * *  in writing.
 * *
 * ********************************************************************************
 * */
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "mkl.h"
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda );
extern void print_int_vector( char* desc, MKL_INT n, MKL_INT* a );
void print_MKL_Complex8( const std::string name, const MKL_Complex8 const * num);
void print_MKL_Complex16( const std::string name, const MKL_Complex16 const * num);
template <class T>
int sumComplex( T * dest, T * src1, T * src2);
template <class T>
int sumComplexArray( T * dest , T * src1, T * one, const size_t nElem);
template <class T>
int addRank(T * a, const size_t nElem_);
template <class T>
void print_MKL_ComplexArray( const std::string name, const T const * arr, const size_t nElem);
/* Parameters */
#define N 4
#define NRHS 2
#define LDA N
#define LDB NRHS

#include <checkpoint.hpp>
#include <cpOptions.h>

/* Main program */
int main(int argc, char* argv[]) {
        MPI_Init(&argc, &argv);
        int myrank, numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        /* Locals */
        MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
        /* Local arrays */
        MKL_INT ipiv[N];
        std::complex<int> c1(1, 2);
        std::complex<double> c2(1.1, 2.1);
//        std::cout << c1 << std::endl << c2 << std::endl;

        MKL_Complex16 one = {0.001, 0.001};
        MKL_Complex8 onef = {0.001f, 0.001f};
//        print_MKL_Complex8("one", &onef);
//        print_MKL_Complex8("aCompf", &aCompf);
        int nElem = 4;
//        MKL_Complex8 a[4] = {
//           { 0.0123f, 0.0550f}, { 0.0791f, 0.0538f}, {0.0980f, 0.0486f}, {0.0732f,  0.057f},
//        };
//        addRank(a, nElem);
        MKL_Complex16 b[4] = {
           {0.033, 0.073}, {0.061, 0.038},
           {0.018, 0.080}, {0.014, 0.071},
        };
        addRank(b, nElem);
        int iteration = 1, ckptIter = 5;

        // ===== Checkpoint definition ===== //
	      Checkpoint myCP( "mklCP", MPI_COMM_WORLD);
        myCP.add("iteration", &iteration);
//        myCP.add("n", &n);
//        myCP.add("c1", &c1);
//        myCP.add("ipiv", ipiv, N);
//        myCP.add("a", a, nElem);
        myCP.add("b", b, nElem);
        myCP.commit(); 
        /* Executable statements */
        printf( "LAPACKE_cgesv (row-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */
        if( myCP.needRestart() ){
		      if(myrank== 0) printf("RESTART ----> \n");
	        myCP.read();
		      iteration++;
	      }
        for (; iteration <= 10; ++iteration)
        {
          MPI_Barrier(MPI_COMM_WORLD);
          if ( myrank == 0) {
            printf("===================== iter: %d ===================== \n", iteration);
          }
//          print_MKL_Complex8Array("a", a, nElem);
          print_MKL_ComplexArray("b", b, nElem);
//          sumComplexArray(a, a , &onef , nElem);
          sumComplexArray(b, b , &one , nElem);
          usleep(500000);
          if(iteration != 0  && (iteration % ckptIter == 0)) {
            printf("====Checkpointing=== \n");
            myCP.update();
            myCP.write();
          }
          

          MPI_Barrier(MPI_COMM_WORLD);
        }
/*
        for (; iteration <= 30; ++iteration){
          MPI_Barrier(MPI_COMM_WORLD);
          sumComplex(&aComp, &aComp, &one);
          sumComplex(&aCompf, &aCompf, &onef);
          MPI_Barrier(MPI_COMM_WORLD);
          usleep(100000);
          if ( myrank == 0) {
            printf("===== iter: %d =====\n", iteration);
          }
          usleep(500000);
          MPI_Barrier(MPI_COMM_WORLD);
//          print_MKL_Complex8("aCompf", &aCompf);
//          print_MKL_Complex16("&&aComp", &aComp);
          computation(a, nElem);
//          info = LAPACKE_cgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
//  Check for the exact singularity
          if( info > 0 ) {
            printf( "The diagonal element of the triangular factor of A,\n" );
            printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
            printf( "the solution could not be computed.\n" );
            exit( 1 );
          }
          if(iteration != 0  && (iteration % ckptIter == 0)) {
            myCP.update();
            myCP.write();
          }

          MPI_Barrier(MPI_COMM_WORLD);
//  Print solution 
//          print_matrix( "Solution", n, nrhs, b, ldb );
//  Print details of LU factorization 
//          print_matrix( "Details of LU factorization", n, n, a, lda );
//  Print pivot indices 
//          print_int_vector( "Pivot indices", n, ipiv );
        }
*/
        MPI_Finalize();
        exit( 0 );
} 

/* End of LAPACKE_cgesv Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.4f,%6.4f)", a[i*lda+j].real, a[i*lda+j].imag );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, MKL_INT n, MKL_INT* a ) {
        MKL_INT j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
        printf( "\n" );
}

void print_MKL_Complex8( const std::string name, const MKL_Complex8 const * num) {
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  printf( "%d_%s: (%6.4f,%6.4f)", myrank, name.c_str(), num->real, num->imag );
  printf( "\n" );
}

void print_MKL_Complex16( const std::string name, const MKL_Complex16 const * num) {
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  printf( "%d_%s: (%6.4f,%6.4f)", myrank, name.c_str(), num->real, num->imag );
  printf( "\n" );
}
template <class T>
int sumComplex( T * dest, T * op1, T * op2){
  dest->real = op1->real + op2->real;
  dest->imag = op1->imag + op2->imag;
  return 0;
}

template <class T>
int sumComplexArray( T * dest , T * src1, T * one, const size_t nElem)
{
  for(int i= 0; i<nElem; ++i)
  {
      dest[i].real = src1[i].real + one->real;
      dest[i].imag = src1[i].imag + one->imag; 
  }
  return 0;
}
template <class T>
int computation(T * a, const size_t nElem){

  return 0;
}

template <class T>
int addRank(T * a, const size_t nElem_){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	for(size_t i= 0; i< nElem_; ++i){
    a[i].real += myrank;
    a[i].imag += myrank;
  }
  return 0;
}

template <class T>
void print_MKL_ComplexArray( const std::string name, const T const * arr, const size_t nElem){
  int myrank, numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Barrier(MPI_COMM_WORLD);
  for( int j = 0; j < numprocs; j++ ) {
    if( j == myrank){
      printf ("%d ===== MKL_ComplexArray: %s =====\n", myrank , name.c_str());
      for( int i = 0; i < nElem; i++ ) {
        std::cout << std::setprecision(6) << std::fixed << "(" << arr[i].real << "," << arr[i].imag << ")" << std::endl;
  
//        printf( " (%6.4f,%6.4f)", arr[i].real, arr[i].imag );
//        printf( "\n" );
      }
    }
    usleep(1000);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return; 
}

