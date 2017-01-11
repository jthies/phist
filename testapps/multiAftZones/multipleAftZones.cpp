// USAGE
// ./minimalWithoutLoop.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40
// In case a failure happens before 40 iterations, the restart flag needs to be provided 
// ./minimalWithoutLoop.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40 -restart 
//

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "multipleAftZones.h"
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

int printflag(int val){
//  if(myrank==0){
    printf("%d:Flag= %d\n", myrank, val);
//  } 
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
	
	Checkpoint cpZ1( "Z1", FT_Comm);
	cpZ1.add("myint", &myint);
	cpZ1.add("iter1", &iter1);
	cpZ1.commit(); 

	if( cpZ1.needRestart()) {
		printf("%d:RESTART Z1 ------>  true \n", myrank);
		cpZ1.read();
    iter1++;
	}
  for(;iter1 <= myCpOpt->getnIter(); ++iter1)
  {
    printflag(myint);
    ++myint;
    usleep(200000); 
    MPI_Barrier(FT_Comm);
    if(iter1 % myCpOpt->getCpFreq() == 0){
		  if(myrank==0){ printf("Checkpointing...\n");}
      cpZ1.update(); 
      cpZ1.write(); 
    }
  }
  
#ifdef AFT
	AFT_END();
#endif

  printf("%d: ============ Z1 END =============== \n", myrank);

// ============ ZONE 2 =============== // 
#ifdef AFT
   	AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 
  CpOptions * myCpOpt = new CpOptions[1];
	MPI_Comm_rank(FT_Comm, &myrank);
  read_params(argc, argv, myCpOpt); 
	MPI_Comm_rank(FT_Comm, &myrank);
	int myint2 = 100;
  int iter2 = 0;
	
	Checkpoint cpZ2( "Z2", FT_Comm);
	cpZ2.add("myint2", &myint2);
	cpZ2.add("iter2", &iter2);
	cpZ2.commit(); 

	if( cpZ2.needRestart() ) {
		printf("%d:RESTART Z2 ------>  true \n", myrank);
		cpZ2.read();
    iter2++;
	}
  for(;iter2 <= myCpOpt->getnIter(); ++iter2)
  {
    printflag(myint2);
    ++myint2;
    usleep(200000); 
    MPI_Barrier(FT_Comm);
    if(iter2 % myCpOpt->getCpFreq() == 0){
		  if(myrank==0){ printf("Checkpointing...\n");}
      cpZ2.update(); 
      cpZ2.write(); 
    }
  }
  
#ifdef AFT
	AFT_END();
#endif






	MPI_Finalize();
	return EXIT_SUCCESS;
}


