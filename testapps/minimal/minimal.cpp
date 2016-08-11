#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "minimal.h"
#include "Checkpoint.hpp"
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
	}
}

int read_params(int argc, char* argv[] , std::string * cpPath){
  prgname = argv[0];

	//=========== Reading commnad line arguments with flags ===============//
	for (int i = 1; i < argc; ++i) {

		if ((!strcmp(argv[i], "-cppath"))) {
			char * tmp = new char[256];
			sprintf(tmp, "%s" ,argv[++i]);
			*cpPath = tmp;
			std::cout << "cpPath: " << *cpPath << std::endl;
		}
	}
	if(cpPath->empty()){
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
   	KSM_MAIN_BEGIN(FT_Comm, &myrank, argv);	
#endif 


	MPI_Comm_rank(FT_Comm, &myrank);
	MPI_Comm_size(FT_Comm, &numprocs);
	printf("%d/%d\n", myrank, numprocs);
	std::string cpPath;
  read_params(argc, argv, &cpPath); 

	int n = 5;
	int myint = 0;
	int * myarray 	= new int[n];
	for(int i = 0; i < n; ++i){
			myarray[i] = 0;
	}

	int iteration = 0, nIter = 20;
	
	Checkpoint * myCP = new Checkpoint[1];
	myCP->setCpPath(cpPath);
	myCP->add("myint", &myint);
	myCP->add("iteration", &iteration);
	myCP->add("myarray", myarray, n);
	myCP->commit(); 
	
	if( failed == true ) {
		failed = false;
		printf("RESTART ------> failed == true \n");
		myCP->read();
	}
  for(; iteration < nIter ; iteration++)
  {
		printf("=== iter: %d\t\n", iteration);
		if(iteration % 5 == 0){
			myCP->update();
			myCP->write();
		}
		myint++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
		usleep(300000);
		MPI_Barrier(FT_Comm);
		if ( iteration+1 == nIter ){
			success = true;
			printf("%d/%d: iterations finishied \n", myrank, numprocs);
	  }
	}
	myCP->update();
	myCP->write();

#ifdef AFT
	KSM_MAIN_END();
#endif

	MPI_Finalize();
	return 0;
}


