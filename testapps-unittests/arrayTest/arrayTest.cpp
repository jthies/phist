#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif


#include "arrayTest.h"
#include <malloc.h>

#include <Checkpoint.hpp>
#include <cp_options.h>

#include <vector>
#include <string>
#include <cstring>
#include <stdio.h>
#include <unistd.h>

static char *prgname = "a.out";

void printusage(){
  int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(myrank==0){		
    printf("Usage: %s [OPTION]...\n",prgname);
    printf("Valid options are:\n");
    printf(" -cppath <PATH-TO-CHECKPOINT>\n");
	}
}

int read_params(int argc, char* argv[] , Cp_Options * myCpOpt){
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
	int myrank, numprocs, printRank = 0;
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

	Cp_Options * myCpOpt = new Cp_Options[1];
  read_params(argc, argv, myCpOpt); 

	int n = 5;
	int * a 	= new int[n];
	double * d 	= new double[n];
	for(int i = 0; i < n; ++i){
			a[i] = 0;
			d[i] = 0.55;
	}
	
	Checkpoint  myCP("CP-L1", myCpOpt->getCpPath(), FT_Comm);
	myCP.add("a", a, n);
	myCP.add("d", d, n);
	myCP.add("iteration", &iteration);
	myCP.commit(); 
	int READ_STATUS = -1;
	if( myCpOpt->getRestartStatus() ) {
		failed = false;
		printf("RESTART ------> failed == true \n");
		READ_STATUS = myCP.read();
		iteration++;
	}
	for(; iteration < myCpOpt->getnIter() ; iteration++)
  {
		printf("=== iter: %d\t a[0]: %d\n", iteration, a[0]);
		for(size_t i = 0; i < n ; ++i){
			a[i] += 1;
			d[i] += 1.0;
		}
		if(iteration % myCpOpt->getCpFreq() == 0){
			myCP.update();
			myCP.write();
		}
		usleep(100000);
		MPI_Barrier(FT_Comm);
		if ( iteration+1 == myCpOpt->getnIter() ){
			success = true;
			printf("%d/%d: iterations finishied \n", myrank, numprocs);
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
		if(myCpOpt->getRestartStatus() == true && iteration == 60 && READ_STATUS == EXIT_SUCCESS)
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST PASSED" << std::endl;	
		}
		else if(myCpOpt->getRestartStatus()==false || READ_STATUS != EXIT_SUCCESS)
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST FAILED: program was not restarted. " << std::endl;	
		}
		else
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST FAILED" << std::endl;	
		}
   } 





#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return 0;
}



