#include "minimal.h"
#include "Checkpoint.hpp"
#include <cstring>
#include <unistd.h>

void read_params(int argc, char* argv[] , std::string * cpPath){

	//=========== Reading commnad line arguments with flags ===============//
	for (int i = 1; i < argc; ++i) {

		if ((!strcmp(argv[i], "-cppath"))) {
			char * tmp = new char[256];
			sprintf(tmp, "%s" ,argv[++i]);
			*cpPath = tmp;
			std::cout << "cpPath: " << *cpPath << std::endl;
		}
	}
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
   	int myrank, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	printf("%d/%d\n", myrank, numprocs);
	std::string cpPath;
  read_params(argc, argv, &cpPath); 
	if(cpPath.empty()){
		printf("cppath not specified.\n");	
	}

	int n = 5;
	int myint = 99;
	int * myarray 	= new int[n];
	for(int i = 0; i < n; ++i){
			myarray[i] = 0;
	}
	
	Checkpoint * myCP = new Checkpoint[1];
	
	myCP->setCpPath(cpPath);
	myCP->add("myint", &myint);
	myCP->add("myarray", myarray, n);
	myCP->commit(); 
 
	int iteration = 0, failed = true, nIter = 20;
	
	int rc = 9;
  for(; iteration < nIter ; iteration++)
  {
		printf("====iter: %d\t\n", iteration);
		if(iteration % 4 == 0){
			myCP->update();
			myCP->write();
		}
		myint++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
		usleep(100000);
	}
   	 
	MPI_Finalize();
	return 0;
}


