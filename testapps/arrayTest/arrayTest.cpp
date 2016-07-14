#include "aft.h"
#include "aft_macros.h"


#include "arrayTest.h"
#include <malloc.h>
#include <Checkpoint.hpp>
#include <vector>
#include <string>
#include <cstring>
#include <stdio.h>

void read_params(int argc, char* argv[] , std::string * cpPath){
	for (int i = 1; i < argc; ++i) {

		if ((!strcmp(argv[i], "-cppath"))) {
			char * tmp = new char[256];
			tmp = argv[++i];
			*cpPath = tmp;
			std::cout << "cpPath: " << *cpPath << std::endl;
		}
	}
}


int main(int argc, char* argv[])
{
	int iteration = 0, failed = true, nIter = 10;
	MPI_Init(&argc, &argv);
   	int myrank, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	printf("%d/%d\n", myrank, numprocs);
	std::string cpPath;
   	read_params(argc, argv, &cpPath); 

	int n = 5;
	int * a 	= new int[n];
	double * d 	= new double[n];
	for(int i = 0; i < n; ++i){
			a[i] = 0;
			d[i] = 0.55;
	}
	
	Checkpoint * myCP = new Checkpoint[1];
	
	myCP->setCpPath(cpPath);
	myCP->add("a", a, n);
	myCP->add("d", d, n);
	myCP->add("iteration", &iteration);
	myCP->commit(); 
 
	printf("%d: START----------------------\n", myrank);
	int rc = 9;
    for(; iteration < nIter ; iteration++)
    {
		printf("===========\n");
		printf("iter: %d\t\n", iteration);
		printf("a[0]: %d\t\n", a[0]);
		for(size_t i = 0; i < n ; ++i){
			a[i] += 1;
			d[i] += 1.0;
		}
/*		if(iteration == 7 && myrank == 1){
			printf("%d_going for kill \n", myrank);	
			exit(0);
		}
*/
		if(iteration % 3 == 0){
			myCP->update();
			myCP->write();
		}
		if(iteration  == 5 && failed == true){
			myCP->read();
			failed = false;
		}
		usleep(1000000);
    }
	return 0;
}



