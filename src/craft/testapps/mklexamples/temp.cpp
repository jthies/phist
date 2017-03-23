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
int sumComplex( T * dest, T * op1, T * op2){
  dest->real = op1->real + op2->real;
  dest->imag = op1->imag + op2->imag;
  return 0;
}
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
        MKL_Complex16 aComp = {0.001, 0.002};
        MKL_Complex8 aCompf = {0.001f, 0.002f};
        aCompf.real += float(myrank);
        aCompf.imag += float(myrank);
        aComp.real += double(myrank);
        aComp.imag += double(myrank);
//        print_MKL_Complex8("one", &onef);
//        print_MKL_Complex8("aCompf", &aCompf);
//        print_MKL_Complex16("aComp", &aComp);
        MKL_Complex8 a[LDA*N] = {
           { 1.23f, -5.50f}, { 7.91f, -5.38f}, {-9.80f, -4.86f}, {-7.32f,  7.57f},
           {-2.14f, -1.12f}, {-9.92f, -0.79f}, {-9.18f, -1.12f}, { 1.37f,  0.43f},
           {-4.30f, -7.10f}, {-6.47f,  2.52f}, {-6.51f, -2.67f}, {-5.86f,  7.38f},
           { 1.27f,  7.29f}, { 8.90f,  6.92f}, {-8.82f,  1.25f}, { 5.41f,  5.37f}
        };
        MKL_Complex8 b[LDB*N] = {
           { 8.33f, -7.32f}, {-6.11f, -3.81f},
           {-6.18f, -4.80f}, { 0.14f, -7.71f},
           {-5.71f, -2.80f}, { 1.41f,  3.40f},
           {-1.60f,  3.08f}, { 8.54f, -4.05f}
        };
        int testint=0, iteration = 1, ckptIter = 5;


        // ===== Checkpoint definition ===== //
	      Checkpoint myCP( "mklCP", MPI_COMM_WORLD);
        myCP.add("testint", &testint);
        myCP.add("iteration", &iteration);
        myCP.add("n", &n);
//        myCP.add("aCompf", &aCompf);
        myCP.add("aComp", &aComp);
//        myCP.add("c1", &c1);
//        myCP.add("ipiv", ipiv, N);
//        myCP.add("a", a, LDA*N);
        myCP.commit(); 
        /* Executable statements */
        printf( "LAPACKE_cgesv (row-major, high-level) Example Program Results\n" );
        /* Solve the equations A*X = B */
        if( myCP.needRestart() ){
		      if(myrank== 0) printf("RESTART ----> \n");
	        myCP.read();
		      iteration++;
	      }
        for (; iteration <= 30; ++iteration){
          MPI_Barrier(MPI_COMM_WORLD);
          sumComplex(&aComp, &aComp, &one);
          MPI_Barrier(MPI_COMM_WORLD);
          usleep(100000);
          if ( myrank == 0) {
            printf("===== iter: %d =====\n", iteration);
          }
          usleep(500000);
          MPI_Barrier(MPI_COMM_WORLD);
//          print_MKL_Complex8("aCompf", &aCompf);
          print_MKL_Complex16("aComp", &aComp);
         
          info = LAPACKE_cgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
          /* Check for the exact singularity */
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
          /* Print solution */
//          print_matrix( "Solution", n, nrhs, b, ldb );
          /* Print details of LU factorization */
//          print_matrix( "Details of LU factorization", n, n, a, lda );
          /* Print pivot indices */
//          print_int_vector( "Pivot indices", n, ipiv );
        }
        MPI_Finalize();
        exit( 0 );
} /* End of LAPACKE_cgesv Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i*lda+j].real, a[i*lda+j].imag );
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
  printf( "%d_%s: (%6.3f,%6.3f)", myrank, name.c_str(), num->real, num->imag );
  printf( "\n" );
}

void print_MKL_Complex16( const std::string name, const MKL_Complex16 const * num) {
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  printf( "%d_%s: (%6.3f,%6.3f)", myrank, name.c_str(), num->real, num->imag );
  printf( "\n" );
}
