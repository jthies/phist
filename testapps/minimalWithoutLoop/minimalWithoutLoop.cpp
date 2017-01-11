// USAGE
// ./minimalWithoutLoop.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40
// In case a failure happens before 40 iterations, the restart flag needs to be provided 
// ./minimalWithoutLoop.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40 -restart 
//

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "minimalWithoutLoop.h"
#include <checkpoint.hpp>
#include <cpOptions.h>
#include <cstring>
#include <unistd.h>


static char *prgname = "a.out";

int myrank=-1;


int read_params(int argc, char* argv[] , CpOptions * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
	std::string cpPathTemp ;
	//=========== Reading commnad line arguments with flags ===============//
	for (int i = 1; i < argc; ++i) {
		if ((!strcmp(argv[i], "-cpfreq"))) {
			sprintf(tmp, "%s" ,argv[++i]);
			myCpOpt->setCpFreq( atoi(tmp) );
			if(myrank==0) std::cout << "cpfreq " << myCpOpt->getCpFreq() << std::endl;
		}
	}
}

int printflag(int val){
  if(myrank==0){
    printf("%d:Flag= %d\n", myrank, val);
  } 
  return 0;
}
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
//	===== AFT BEGIN =====
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
   	AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 
	CpOptions * myCpOpt = new CpOptions[1];
	MPI_Comm_rank(FT_Comm, &myrank);
  read_params(argc, argv, myCpOpt); 

	int n = 5;
	int myint = 0;
  int iter1 = 0;
	
	Checkpoint myCP( "minimalWithoutloop", FT_Comm);
	myCP.add("myint", &myint);
	myCP.commit(); 
	
	if( myCP.needRestart()) 
  {
		printf("%d:RESTART ------>  true \n", myrank);
		myCP.read();
	}
  sleep(5);
  MPI_Barrier(FT_Comm);
  printf("%d_F1 \n", myrank);
  
  sleep(5);

  MPI_Barrier(FT_Comm);
  printf("%d_F2 \n", myrank);
  sleep(5);
  MPI_Barrier(FT_Comm);
  printflag(3);
 /* 
  for(; iter1 <= myCpOpt->getnIter() ; iter1++)
  {
		myint++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
		usleep(300000);
		{ printf("=== iter: %d , myint: %d \t\n", iter1, myint-1);}
		if(iter1 % myCpOpt->getCpFreq() == 0){
		  if(myrank==0){ printf("Checkpointing...\n");}
			myCP.update();
			myCP.write();
		}
		MPI_Barrier(FT_Comm);
	}
*/
#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return EXIT_SUCCESS;
}


