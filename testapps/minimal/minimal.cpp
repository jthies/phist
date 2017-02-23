// USAGE
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40
// In case a failure happens before 40 iterations, the restart flag needs to be provided 
// ./minimal.bin -cppath <CHECKPOINT-PATH> -cpfreq 10 -niter 40 -restart 
//

#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "minimal.h"
#include <craft.h>
#include <cstring>
#include <unistd.h>

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

static char *prgname = "a.out";

int read_params(int argc, char* argv[] , CpOptions * myCpOpt){
  prgname = argv[0];
	char * tmp = new char[256];
  int myrank=-1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
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
  int myrank, numprocs;
//	===== AFT BEGIN =====
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
  AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 
	CpOptions * myCpOpt = new CpOptions[1];
  read_params(argc, argv, myCpOpt); 

	int n = 5;
	int myint = 0;
  long double ld = 0.99;
  std::complex<double> d_c1(0.1, 0.1);
  std::complex<double> d_one(1.0, 1.0);
  std::complex<int> i_c1(0, 0);
  std::complex<int> i_one(1, 1);
//  using namespace std::complex_literals;
//  std::cout << std::fixed << std::setprecision(1);
//  std::complex<double> z1 = 1i * 1i;     // imaginary unit squared
//  std::cout << "i * i = " << z1 << '\n'; 
	double mydouble = 0.0123;
	int * myarray 	= new int[n];
	for(int i = 0; i < n; ++i){
			myarray[i] = 0;
	}

	int iteration = 0;
	
	Checkpoint myCP( "aa", FT_Comm);
	myCP.add("myint", &myint);
//	myCP.add("mydouble", &mydouble);
	myCP.add("iteration", &iteration);
  myCP.add("myarray", myarray, n);
//  myCP.add("ld", &ld);
//  myCP.add("i_c1", &i_c1);
//  myCP.add("d_c1", &d_c1);
	myCP.commit(); 
	
  if( myCP.needRestart() == true) 
  {
    { if (myrank == 0) printf("%d:RESTART ------>  true \n", myrank); }
    if(myCP.read()==EXIT_SUCCESS){
		  iteration++;
		  { printf("%d:iteration = %d \n", myrank, iteration); }
    }
	}
  for(; iteration <= myCpOpt->getnIter() ; iteration++)
  {
		myint++;
		mydouble++;
		for(size_t i = 0; i < n ; ++i){
			myarray[i] += 1;
		}
//    if ( iteration == 15) {
//      char * cmd=new char[256];
//      sprintf(cmd, "$HOME/killproc.sh 2 mini");
//      if (myrank == 0) {
//        system(cmd); 
//      }
//    }
		{ printf("%d:=== iter: %d , myint: %d \t\n", myrank, iteration, myint-1);}
//		{ printf("%d:=== iter: %d , c1: %d \t\n", myrank, iteration, c1);}
		usleep(1000000);
//		{ if(myrank==0) printf("%d:=== iter: %d , mydouble: %f \t\n", myrank, iteration, mydouble-1);}
		if(iteration % myCpOpt->getCpFreq() == 0){
		  if(myrank==0){ printf("Checkpointing...\n");}
			myCP.update();
  		myCP.write();
		}
	  MPI_Barrier(FT_Comm);
	}
#ifdef AFT
	AFT_END()
#endif

	MPI_Finalize();
	return EXIT_SUCCESS;
}


