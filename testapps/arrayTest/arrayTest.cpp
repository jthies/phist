#include <vector>
#include <string>
#include <cstring>
#include <stdio.h>
#include <unistd.h>
#include <malloc.h>

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif
#include <checkpoint.hpp>
#include <cpOptions.h>

#include "arrayTest.h"

static char *prgname = "a.out";

template < typename T>
void printArray( T * ptr, const int nRows){
	for(size_t i = 0; i < nRows; ++i){
		std::cout << ptr[i] << std::endl;
	}
	return;
}

template < typename T>
void printDoubleArray( T ** ptr, const int nRows, const int nCols){
	for(size_t i = 0; i < nRows; ++i){
		for(size_t j = 0; j < nCols; ++j){
			std::cout << ptr[j][i] << "\t";
		}
		std::cout << std::endl;
	}
	return;
}
void printusage(){
  int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(myrank==0){		
    printf("Usage: %s [OPTION]...\n",prgname);
    printf("Valid options are:\n");
    printf(" -cppath <PATH-TO-CHECKPOINT>\n");
	}
}

int read_params(int argc, char* argv[] , CpOptions * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
	std::string cpPathTemp ;
	for (int i = 1; i < argc; ++i) {
		if ((!strcmp(argv[i], "-cppath"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setCpPath(tmp);
			cpPathTemp = myCpOpt->getCpPath();
			std::cout << "cpPath: " << cpPathTemp << std::endl;
		}
		if ((!strcmp(argv[i], "-restart"))) {
			bool restart = true;
			myCpOpt->setRestartStatus( restart );
			std::cout << "Restart " << restart << std::endl;
		}
		if ((!strcmp(argv[i], "-niter"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setnIter( atoi(tmp) );
			std::cout << "nIter " << myCpOpt->getnIter() << std::endl;
		}
		if ((!strcmp(argv[i], "-cpfreq"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setCpFreq( atoi(tmp) );
			std::cout << "cpfreq " << myCpOpt->getCpFreq() << std::endl;
		}
	}
	if(cpPathTemp.empty()){
		printusage();
		exit(EXIT_FAILURE);
	}
}


int main(int argc, char* argv[])
{
	int iteration = 0;
	MPI_Init(&argc, &argv);
  int myrank, numprocs;
//	===== AFT BEGIN =====
	int success = false;
	int failed = false;
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
   	AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 

	MPI_Comm_rank(FT_Comm, &myrank);
	MPI_Comm_size(FT_Comm, &numprocs);
	printf("%d/%d\n", myrank, numprocs);

	CpOptions * myCpOpt = new CpOptions[1];
  read_params(argc, argv, myCpOpt); 

	int n = 5, nRows=5, nCols=5;
	printf("sizeof(float)=%d\n", sizeof(float));

	int * a 	= new int[n];
	float ** aa = new float*[nCols];
	for(int i=0; i<nRows; ++i){
		aa[i] = new float[nRows];	
	} 
	for(int i=0; i<nCols; ++i){
		for(int j=0; j<nRows; ++j){
			aa[i][j] = j+(0.1*i);
		}
	}
	printDoubleArray(aa, nRows, nCols);
	double * d 	= new double[n];
	for(int i = 0; i < n; ++i){
			a[i] = 0;
			d[i] = 0.55;
	}
	int rankCP = myrank;
	Checkpoint  myCP( myCpOpt->getCpPath(), FT_Comm);
  
	myCP.add("a", a, n);
	myCP.add("aa", aa, nRows, nCols, ALL);
	myCP.add("d", d, n);
	myCP.add("iteration", &iteration);
	myCP.add("rankCP", &rankCP);
	myCP.commit(); 

	if( myCpOpt->getRestartStatus() ) {
		printf("RESTART ------> failed == true \n");
		failed = false;
		myCP.read();
		iteration++;
	}

	for(; iteration <= myCpOpt->getnIter() ; iteration++)
  {
		printf("=== iter: %d\t a[0]: %d\t aa[0][0]:%f\n", iteration, a[0], aa[0][0]);
		for(size_t i = 0; i < n ; ++i){
			a[i] += 1;
			d[i] += 1.0;
		}
		for(int i=0; i<nCols; ++i){
			for(int j=0; j<nRows; ++j){
				aa[i][j] += 0.001;
			}
		}
/*		if(iteration == 7 && myrank == 1){
			printf("%d_going for kill \n", myrank);	
			exit(0);
		}
*/
		if(iteration % myCpOpt->getCpFreq() == 0){
			myCP.update();
			myCP.write();
		}
		usleep(200000);
		MPI_Barrier(FT_Comm);
		if ( iteration+1 == myCpOpt->getnIter() ){
			success = true;
			printf("%d/%d: iterations finishied \n", myrank, numprocs);
	  }
  }
//	printf("rankCP: %d\n", rankCP);
	printArray(a, n);
	printArray(d, n);
	printDoubleArray(aa , nRows, nCols);
#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return 0;
}



