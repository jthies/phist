// USAGE
// mpirun -am ft-enable-mpi --mca btl openib,sm,self -np 2 -npernode 1 ./minimalWithoutLoop.bin
//

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include <craft.h>
#include <cstring>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#include "minimalWithoutLoop.h"

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
	int n = 5;
  int restartedRun = false;
	MPI_Init(&argc, &argv);
//	===== AFT BEGIN =====
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
   	AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 
	MPI_Comm_rank(FT_Comm, &myrank);
  Checkpoint CP1("miniWLoop", FT_Comm);
  CP1.add("n", &n);
  CP1.commit(); 

  if(CP1.needRestart()){
    restartedRun = true;
  }

  if(restartedRun == false){
    CP1.update(); CP1.write();
    if(myrank==0) printf("=== 1st run:Start ===\n");
    MPI_Barrier(FT_Comm);
    for(int i=0; i<10000000; i++) {
      double res=0.0, res1=0.0; 
      if ( i == 15) {
        char * cmd=new char[256];
        sprintf(cmd, "$HOME/killproc.sh 2 mini");
        if (myrank == 0) {
          system(cmd); 
        }
      }
      for(int j=0;j<10000000;j++){
        int rand1= rand();
        int rand2= rand();
        res += (double)rand1/(double)rand2 + 2.00;
      }
      if (res < 2.00) break; 
      if(myrank==0) printf("i:%d _ %f\n", i, res);
      MPI_Barrier(FT_Comm);
    }
  }
  if(restartedRun == true ){
    if(myrank==0) printf("=== 2st run ===\n");
    MPI_Barrier(FT_Comm);
    sleep(5);
    MPI_Barrier(FT_Comm);
  }
  
#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return EXIT_SUCCESS;
}


