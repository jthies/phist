#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif


#include "arrayTest.h"
#include <malloc.h>

#include <checkpoint.hpp>
#include <cpOptions.h>

#include <vector>
#include <string>
#include <cstring>
#include <stdio.h>
#include <unistd.h>

static char *prgname = "a.out";

int read_params(int argc, char* argv[] , CpOptions * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
	for (int i = 1; i < argc; ++i) {
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
}

int main(int argc, char* argv[])
{
	int iteration = 0;
	MPI_Init(&argc, &argv);
	int myrank, numprocs, printRank = 0;
//	===== AFT BEGIN =====
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

	int n = 5;
	int * a 	= new int[n];
	double * d 	= new double[n];
	for(int i = 0; i < n; ++i){
			a[i] = 0;
			d[i] = 0.55;
	}
	
	Checkpoint  myCP("arraytest", FT_Comm);
	myCP.add("a", a, n);
	myCP.add("d", d, n);
	myCP.add("iteration", &iteration);
	myCP.commit(); 
	int READ_STATUS = -1;
	if( myCP.needRestart()) {
		printf("RESTART ------> \n");
		READ_STATUS = myCP.read();
		iteration++;
	}
	for(; iteration <= myCpOpt->getnIter() ; iteration++)
  {
		printf("=== iter: %d\t a[0]: %d\n", iteration, a[0]);
		for(size_t i = 0; i < n ; ++i){
			a[i] += 1;
			d[i] += 1.0;
		}
		usleep(100000);
		MPI_Barrier(FT_Comm);
		if(iteration % myCpOpt->getCpFreq() == 0){
			myCP.update();
			myCP.write();
		}
  }
	for(int i = 0; i < n; ++i){
		if(myrank==printRank) 
			printf("a[%d]: %d \n", i, a[i]);
	}
	for(int i = 0; i < n; ++i){
		if(myrank==printRank) 
			printf("d[%d]: %f \n", i, d[i]);
	}

	MPI_Barrier(FT_Comm);
	{
	  if(myrank==printRank) 
   	printf("-------------------------------------\n");
		if(iteration-1 == 60 && READ_STATUS == EXIT_SUCCESS)
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST PASSED" << std::endl;	
		}
		else if(iteration-1 != 60 && READ_STATUS != EXIT_SUCCESS)
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST FAILED: program was not restarted or READ_STATUS==EXIT_FAILURE" << std::endl;	
		}
		else
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST FAILED" << "iteration: " << iteration << "READ_STATUS: " << READ_STATUS << std::endl;	
		}
  } 


#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return 0;
}



