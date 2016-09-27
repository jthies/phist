// USAGE
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40
// In case a failure happens before 40 iterations, the restart flag needs to be provided 
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40 -restart 
//

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "minimal.h"
#include <Checkpoint.hpp>
#include <cp_options.h>
#include <cstring>
#include <unistd.h>


static char *prgname = "a.out";

void printusage(){
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if(myrank==0){	
    printf("Usage: %s [OPTION]...\n",prgname);
    printf("Valid options are:\n");
    printf(" -cppath <PATH-TO-CHECKPOINT>\n");
    printf(" -restart : In case of a restart\n");
	}
}

int read_params(int argc, char* argv[] , Cp_Options * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
	std::string cpPathTemp ;
	//=========== Reading commnad line arguments with flags ===============//
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

	Cp_Options * myCpOpt = new Cp_Options[1];
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
	
	Checkpoint  myCP("CP-L1", myCpOpt->getCpPath(), FT_Comm);
//	myCP.disableSCR();
	myCP.add("myint", &myint);
	myCP.add("mydouble", &mydouble);
	myCP.add("iteration", &iteration);
	myCP.add("myarray", myarray, n);
	myCP.commit(); 
	
	if( myCpOpt->getRestartStatus() ) {
		failed = false;
		printf("RESTART ------> failed == true \n");
		myCP.read();
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
		printf("=== iter: %d , myint: %d \t\n", iteration, myint-1);
		if(iteration % myCpOpt->getCpFreq() == 0){
			myCP.update();
			myCP.write();
		}
		
		MPI_Barrier(FT_Comm);
		if ( iteration+1 == myCpOpt->getnIter() ){
			success = true;
			printf("%d/%d: iterations finishied \n", myrank, numprocs);
	  }
	}

#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return 0;
}


