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
#include <checkpoint.hpp>
#include <cpOptions.h>
#include <cstring>
#include <unistd.h>


static char *prgname = "a.out";


int read_params(int argc, char* argv[] , CpOptions * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
	std::string cpPathTemp ;
  int myrank=-1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//=========== Reading commnad line arguments with flags ===============//
	for (int i = 1; i < argc; ++i) {
		if ((!strcmp(argv[i], "-niter"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setnIter( atoi(tmp) );
			if(myrank==0) std::cout << "nIter " << myCpOpt->getnIter() << std::endl;
		}
		if ((!strcmp(argv[i], "-cpfreq"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setCpFreq( atoi(tmp) );
			if(myrank==0) std::cout << "cpfreq " << myCpOpt->getCpFreq() << std::endl;
		}
	}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
  int myrankTest, numprocs;
//	===== AFT BEGIN =====
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
  AFT_BEGIN(FT_Comm, &myrankTest, argv);	
#endif 
	CpOptions * myCpOpt = new CpOptions[1];
  read_params(argc, argv, myCpOpt); 

	int n = 5;
	int myint = 0;
	double mydouble = 0.0123;
	int * myarray 	= new int[n];
	for(int i = 0; i < n; ++i){
			myarray[i] = 0;
	}

	int iteration = 0;
	
	Checkpoint myCP( "a", FT_Comm);
	myCP.disableSCR();
	myCP.add("myint", &myint);
	myCP.add("mydouble", &mydouble);
	myCP.add("iteration", &iteration);
	myCP.add("myarray", myarray, n);
	myCP.commit(); 
	
//	if( myCP.getRestartStatus() == true) 
  if( myCP.needRestart() == true) 
  {
    printf("%d:RESTART ------>  true \n", myrankTest);
    if(myCP.read()==EXIT_SUCCESS){
		  iteration++;
		  printf("%d:iteration = %d \n", myrankTest, iteration);
    }
	}
  for(; iteration <= myCpOpt->getnIter() ; iteration++)
  {
		myint++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
		usleep(300000);
		{ printf("=== iter: %d , myint: %d \t\n", iteration, myint-1);}
		if(iteration % myCpOpt->getCpFreq() == 0){
		  if(myrankTest==0){ printf("Checkpointing...\n");}
			myCP.update();
			myCP.write();
		}
	  MPI_Barrier(FT_Comm);
	}
#ifdef AFT
	AFT_END()
#endif

	MPI_Finalize();
	return EXIT_SUCCESS;
}


