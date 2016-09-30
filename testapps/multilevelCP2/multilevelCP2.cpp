// USAGE
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40
// In case a failure happens before 40 iterations, the restart flag needs to be provided 
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40 -restart 
//

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "multilevelCP2.h"
#include <Checkpoint.hpp>
#include <cstring>
#include <unistd.h>


static char *prgname = "a.out";

template < typename T>
void printArray( T * ptr, const int nRows){
	for(size_t i = 0; i < nRows; ++i){
		std::cout << ptr[i] << std::endl;
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
    printf(" -restart : In case of a restart\n");
	}
}

int read_params(int argc, char* argv[] , Cp_Options * myCpOpt){
  prgname = argv[0];
	printf("====== argv[0] ====== %s\n", argv[0]);
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

	int n = 2;
	double a 	= 0;
	double * b 	= new double[n];
	for(int i = 0; i < n; ++i){
			b[i] = 0.0;
	}
	
	int naIter = myCpOpt->getnIter();
	int aIter=0;
	
//	Checkpoint  cpL1;
	Checkpoint  cpL1(myCpOpt->getCpPath(), FT_Comm);
	cpL1.add("a", &a);
	cpL1.add("aIter", &aIter);
	cpL1.commit(); 

	if( myCpOpt->getRestartStatus() ) {
		printf("RESTART L1 ------> failed == true \n");
		cpL1.read();
		aIter++;
		printf("aIter to RESTART = %d \n", aIter);
	}
  for(; aIter < naIter; aIter++)
  {
			int nbIter = 10, bIter = 0;
			static Checkpoint  cpL2(myCpOpt->getCpPath(), FT_Comm);
			cpL2.add("b", b);
			cpL2.add("bIter", &bIter);
			cpL2.commit(); 
	
			bIter = 0;
			if( myCpOpt->getRestartStatus() ) {
				bool restart = false;
				myCpOpt->setRestartStatus( restart );
				printf("RESTART L2------> failed == true \n");
				cpL2.read();
				bIter++;
				printf("bIter to RESTART = %d \n", bIter);
			}
  		for(double x=0.0; bIter < nbIter ; bIter++,x++)
  		{
				for(size_t i = 0; i < n ; ++i){
					b[i] += 0.001;
				}
				usleep(100000);
				printf("=== aIter: %d ,bIter: %d, \tb[0]: %f \n", aIter, bIter, b[0] );
//				if(bIter % 5 == 0){
				if(bIter % 5 == 0 ){						// if inner checkpoint is not saved on the last iteration and new data depends on the previous data. data after reading will be corrupted.
					cpL2.update();
					cpL2.write();
				}
			}
			a += b[0];

			usleep(100000);
			printf("================================================ aIter: %d, \ta[0]: %f \n", aIter, a );
			if(nbIter % 1 == 0){
					cpL1.update();
					cpL1.write();
			}
	}

#ifdef AFT
	AFT_END();
#endif

	MPI_Finalize();
	return 0;
}


