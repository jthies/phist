// USAGE
// ./minimal.bin -cpfreq 10 -niter 40

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "minimal.h"
#include <checkpoint.hpp>
#include <cpOptions.h>
#include <cstring>
#include <unistd.h>

int printRank = 0;
static char *prgname = "a.out";
int myrank;
int read_params(int argc, char* argv[] , CpOptions * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
	//=========== Reading commnad line arguments with flags ===============//
	for (int i = 1; i < argc; ++i) {
		if ((!strcmp(argv[i], "-niter"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setnIter( atoi(tmp) );
			if(myrank==printRank) std::cout << "nIter " << myCpOpt->getnIter() << std::endl;
		}
		if ((!strcmp(argv[i], "-cpfreq"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setCpFreq( atoi(tmp) );
			if(myrank==printRank) std::cout << "cpfreq " << myCpOpt->getCpFreq() << std::endl;
		}
	}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
  int numprocs;

//	===== AFT BEGIN =====
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
   	AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 

	CpOptions * myCpOpt = new CpOptions[1];
	MPI_Comm_rank(FT_Comm, &myrank);
	MPI_Comm_size(FT_Comm, &numprocs);
  read_params(argc, argv, myCpOpt); 

	int n = 5;
	int myint = 0;
	double mydouble = 0.0123;
	int * myarray 	= new int[n];
	for(int i = 0; i < n; ++i){
			myarray[i] = 0;
	}

	int iteration = 0;
	
	Checkpoint  myCP("mini", FT_Comm);
//	myCP.disableSCR();
	myCP.add("myint", &myint);
	myCP.add("iteration", &iteration);
	myCP.add("myarray", myarray, n);
	myCP.commit(); 
	

	int READ_STATUS=-1;
	if( myCP.needRestart() ) {
		printf("RESTART ---------> \n");
		READ_STATUS = myCP.read();
		iteration++;
		printf("iteration = %d \n", iteration);
	}
  for(; iteration <= myCpOpt->getnIter() ; iteration++)
  {
		myint++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
		usleep(100000);
		MPI_Barrier(FT_Comm);
	  if(myrank==printRank) 
		  printf("=== iter: %d , myint: %d \t\n", iteration, myint-1);
		if(iteration % myCpOpt->getCpFreq() == 0){
			myCP.update();
			myCP.write();
		}
	}
	{
	  if(myrank==printRank) 
   	printf("-------------------------------------\n");
		if( iteration-1 == 60 && READ_STATUS == EXIT_SUCCESS )
		{
			if(myrank==printRank) 
				std::cout << "CRAFT TEST PASSED" << std::endl;	
		}
		else if( iteration-1 !=60 && READ_STATUS != EXIT_SUCCESS)
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


