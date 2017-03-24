// USAGE
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40
// In case a failure happens before 40 iterations, the restart flag needs to be provided 
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40 -restart 
//
#include <cstring>
#include <unistd.h>

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "multilevelCP.h"
#include <checkpoint.hpp>

static char *prgname = "a.out";
int myrank;

template < typename T>
void printArray( T * ptr, const int nRows){
	for(size_t i = 0; i < nRows; ++i){
		if(myrank==0) std::cout << ptr[i] << std::endl;
	}
	return;
}


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

	int n = 2;
	double a 	= 0.;
	double * b 	= new double[n];
	for(int i = 0; i < n; ++i){
			b[i] = 0.0;
	}
	
	int naIter = myCpOpt->getnIter();
	int nbIter = 50;
	int aIter=0, bIter=0;
	
	Checkpoint  cpL1("cpL1", FT_Comm);
	Checkpoint  cpL2("cpL2", FT_Comm);
  cpL2.disableSCR();
  cpL1.disableSCR();
	cpL1.add("a", &a);
	cpL1.add("aIter", &aIter);
	cpL1.commit(); 
	cpL2.add("b", b);
	cpL2.add("bIter", &bIter);
	cpL2.commit(); 
	if( cpL1.needRestart() == true) {
	  printf("RESTART L1------> failed == true \n");
    if(cpL1.read()==EXIT_SUCCESS){
		  aIter++;
		  if(myrank==0) std::cout << "aIter =" << aIter << std::endl;
    }
	}
  printf("%d: before for loop\n", myrank );
  for(; aIter < myCpOpt->getnIter(); aIter++)
  {
			std::cout << "aIter:" << aIter << std::endl;
			bIter = 0;
			if( cpL2.needRestart() == true ) {
				printf("RESTART L2------> failed == true \n");
				if(cpL2.read()==EXIT_SUCCESS){
				  bIter++;
				  printf("bIter = %d \n", bIter);
        }
			}
  		for(; bIter < nbIter ; bIter++)
  		{
				for(size_t i = 0; i < n ; ++i){
					b[i] += 0.1;
				}
				usleep(200000);
        MPI_Barrier(FT_Comm);
        printf("=== bIter: %d, \tb[0]: %f \n", bIter, b[0] );
				if( (bIter % myCpOpt->getCpFreq()) == 0){
          printf("Checkpointing b...\n");
					cpL2.update();
					cpL2.write();
				}
			}
			a = b[0];
			usleep(200000);
      printf("============ aIter: %d, \ta: %f \n", aIter, a);
			if(aIter % 1 == 0){
          printf("Checkpointing a...\n");
					cpL1.update();
					cpL1.write();
			}

	}
  if(myrank==0) std::cout <<("1before AFT_END\n");
#ifdef AFT
	AFT_END();
#endif

  if(myrank==0) std::cout <<"before MPI_Finalize\n";
	MPI_Finalize();
	return 0;
}


