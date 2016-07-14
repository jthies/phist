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

	int i1 = 99;
	double d1 = 88.5; 
	int n = 5;
	int * a 	= new int[n];
	double * d 	= new double[n];
	for(int i = 0; i < n; ++i){
			a[i] = 0;
			d[i] = 0;
	}
	
	Checkpoint * myCP = new Checkpoint[1];
	
	myCP->setCpPath(cpPath);
	myCP->add("i1", &i1);
	myCP->add("d1", &d1);
//	CpPOD<int>  * myCpPOD1= new CpPOD<int>[1];	
//	CpPOD<double>  * myCpPOD2= new CpPOD<double>[1];	
//	myCpPOD1->init(&i1);
//	myCpPOD2->init(&d1);
	
//	myCP->add("myCpPOD1", myCpPOD1);
//	myCP->add("myCpPOD2", myCpPOD2);
//	myCP->commit(); 
 
	printf("%d: START----------------------\n", myrank);
	int iteration = 0, failed = true, nIter = 10;
/*	if(failed == true){
		failed = false;
		printf("RESTART ----> failed == true \n");
		myCP->read();
		iteration++;
	
	}
	*/
	int rc = 9;
    for(; iteration < nIter ; iteration++)
    {
		printf("===========\n");
		printf("iter: %d\t\n", iteration);
//		myA1.print();
//		printf("a[0]: %d\t\n", a[0]);
/*		if(iteration == 7 && myrank == 1){
			printf("%d_going for kill \n", myrank);	
			exit(0);
		}
*/
		//rc = MPI_Barrier(MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD);
//		printf("iteration: %d\n");
		if(iteration % 3 == 0){
			myCP->update();
			myCP->write();
		}
		if(iteration  == 5){
	//		myCP->read();
		}
		i1++;
		d1 = d1 - 1.0;
//		myA1.print();
		usleep(1000000);
    }
   	 
    printf("%d: ------------------------------------- END\n", myrank);
	MPI_Finalize();
	return 0;
}

/*  //UserClass *my_object = new UserClass();
  int my_int=42;
  double my_double=42.9;
  Checkpoint backup;
  backup.add("i",&my_int);
  backup.add("x",&my_double);
 // backup.add("my_object",my_object);
  backup.write();
}
*/

