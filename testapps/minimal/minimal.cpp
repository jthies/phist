#include "minimal.h"
#include "Checkpoint.hpp"

//#include "UserClass.hpp"

void read_params(int argc, char* argv[] , std::string * cpPath){

	//============================ Reading commnad line arguments with flags ===============//
	for (int i = 1; i < argc; ++i) {

		if ((!strcmp(argv[i], "-cppath"))) {
			char * tmp = new char[256];
			tmp = argv[++i];
			*cpPath = tmp;
			std::cout << "cpPath: " << *cpPath << std::endl;
		}
/*		if ((!strcmp(argv[i], "-m"))) {
			m = atoi(argv[++i]);
			printf("m: %d\n",m);
		}
*/
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
 
	printf("%d: START----------------------\n", myrank);
	int iteration = 0, failed = true, nIter = 10;
	
	int rc = 9;
    for(; iteration < nIter ; iteration++)
    {
		printf("===========\n");
		printf("iter: %d\t\n", iteration);
		if(iteration % 3 == 0){
			myCP->update();
			myCP->write();
		}
		if(iteration  == 5){
	//		myCP->read();
		}
		myint++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
		usleep(1000000);
    }
   	 
    printf("%d: ------------------------------------- END\n", myrank);
	MPI_Finalize();
	return 0;
}


