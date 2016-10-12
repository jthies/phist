#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

#include "ghostTest.h"
#include <malloc.h>
#include <checkpoint.hpp>
#include <cp_options.h>

#include <vector>
#include <string>
#include <cstring>
#include <typeinfo>
#include <ghost.h>
#include <ghost/types.h>
#include <stdio.h>
extern "C"{
#include "essexamples.h"
}
#include <math.h>

#include <cpGhost.hpp>
#include <cpPOD.hpp>

#define DP
#ifdef DP
    GHOST_REGISTER_DT_D(vecdt) // vectors have double values
    GHOST_REGISTER_DT_D(matdt) // matrix has double values
#define EPS 1e-32
#else
    GHOST_REGISTER_DT_S(vecdt) // vectors have float values
    GHOST_REGISTER_DT_S(matdt) // matrix has float values
#define EPS 1e-16
#endif

#define N 4

typedef struct args{
	int argc;
	char **argv;
}args;

static int diag(ghost_gidx row, ghost_lidx *rowlen, ghost_gidx *col, void *val, __attribute__((unused)) void *arg)
{ 
	*rowlen = 1; 
	col[0] = row; 
	//((mydata_t *)val)[0] = (mydata_t)(row+1); 
	((double *)val)[0] = 2.0;
	return 0; 
}

static void *mainTask(void *varg)
{
	args * arg = (args * ) varg;
	int argc = arg->argc;
	char **argv = arg->argv;
	
	int printRank=0;
  int myrank, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	printf("%d/%d\n", myrank, numprocs);

	int success = false;
	int failed = false;
  MPI_Comm FT_Comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FT_Comm);
#ifdef AFT
  AFT_BEGIN(FT_Comm, &myrank, argv);	
#endif 

    ghost_context *context;
	
    matdt_t alpha = 0., alphaold = 0., lambda=0., lambdaold = 0., zero = 0., neg = -1., one = 1., tmp = 0.;
    int iteration = 0, nIter = 1000, cp_freq = 99999, restart = -1;
	char * cpPath = new char[256];
    double start = 0.;
    matdt_t localdot[3];

    ghost_sparsemat *A;
    ghost_densemat *x;
    ghost_densemat *r;
    ghost_densemat *b;
    ghost_densemat *p;
    ghost_densemat *v;
    ghost_densemat **vv;
	
    essexamples_get_iterations(&nIter);
		essexamples_get_cp_freq(&cp_freq);
		essexamples_get_cp_folder(&cpPath);
		essexamples_get_restart(&restart);

    ghost_densemat_traits vtraits = GHOST_DENSEMAT_TRAITS_INITIALIZER;
    vtraits.datatype = vecdt;
    	
    ghost_sparsemat_traits mtraits = GHOST_SPARSEMAT_TRAITS_INITIALIZER;
#ifdef AFT
    essexamples_create_context_and_matrix_ft(&context,&A,&mtraits, FT_Comm);
#else 
    essexamples_create_context_and_matrix(&context,&A,&mtraits);
#endif
    ghost_spmv_flags spmvmOpt = GHOST_SPMV_DOT;
    ghost_spmv_opts spmvtraits = GHOST_SPMV_OPTS_INITIALIZER;
    spmvtraits.dot = localdot;
    spmvtraits.flags = GHOST_SPMV_DOT;

		ghost_rank(&myrank, context->mpicomm);

    essexamples_create_densemat(&x,&vtraits,context);
    essexamples_create_densemat(&b,&vtraits,context);
    essexamples_create_densemat(&r,&vtraits,context);
    essexamples_create_densemat(&v,&vtraits,context);
    essexamples_create_densemat(&p,&vtraits,context);
	
    ghost_densemat_init_val(x,&one);  // x = 1, start with defined value
    ghost_densemat_init_val(r,&zero); // r = 0
    ghost_densemat_init_val(b,&zero); // b = 0
    ghost_densemat_init_val(v,&zero); // v = 0

    ghost_spmv(r,A,x,spmvtraits); // r = A * x
    ghost_axpby(r,b,&one, &neg);  // r = -1*r + b
    ghost_dot(&alpha,r,r,A->context->mpicomm);        // alpha = ||r||

    ghost_densemat_init_densemat(p,r,0,0); // p = r

		if( myrank == printRank) {
				printf("==== Defining CP ====\n");
		}
	Checkpoint  myCP( cpPath, FT_Comm);
	myCP.disableSCR();
	myCP.add("iteration", &iteration);	
	myCP.add("lambda", &lambda);	
	myCP.add("alpha", &alpha);	
	myCP.add("x", x);
	myCP.add("r", r);
//	myCP.add("vv", vv, 4);
	myCP.add("p", p);
	// this class object contains only an integer
  myCP.commit();
    	
	iteration = 0;
	if(restart == true){
		failed = false;
		if(myrank== printRank) printf("RESTART ----> failed == true \n");
		myCP.read();
		iteration++;
	}
	
    for(; iteration < nIter && alpha > EPS; iteration++)
    {
        alphaold = alpha;
        lambdaold = lambda;

        localdot[0] = 0.;
        localdot[1] = 0.;
        localdot[2] = 0.;
        ghost_spmv(v,A,p,spmvtraits); // v = A * p

        lambda = localdot[1];

        lambda = alphaold / lambda; // ...

        ghost_instr_prefix_set("x_");
        ghost_axpy(x,p,&lambda); // x = x + lambda*p

        ghost_instr_prefix_set("r_");
        lambda *= -1.;
        ghost_axpy(r,v,&lambda); // r = r + lambda*v

        ghost_instr_prefix_set("");
        lambda *= -1.;

        ghost_dot(&alpha,r,r,A->context->mpicomm); // alpha = ||r||

        tmp = alpha/alphaold;

        ghost_axpby(p,r,&one,&tmp);			//	p = 1*r + tmp*p // ghost_axpby(y, x, a, b)
		
				usleep(1000);				// TODO: just for testing, in order to see the progress of program
        if (myrank == printRank){
					printf("iter=%d, ", iteration);
          printf("alpha[%4d] = %g\n",iteration+1, alpha);
	   		}
				if(iteration % cp_freq == 0){
					if(myrank == printRank ) 
							printf("iteration=%d, cp_freq=%d\n", iteration, cp_freq);
					myCP.update();
					myCP.write();
				}
        		//fflush(stdout);
	   		if ( iteration+1 == nIter || alpha <= EPS){
					success = true;
					printf("%d/%d: iterations finishied \n", myrank, numprocs);
	   		}
    }
    
    // print norm of residual
    r->fromScalar(r,&zero);
    ghost_spmv(r,A,x,GHOST_SPMV_OPTS_INITIALIZER);
    r->axpy(r,b,&neg);
    vecdt_t rnorm;
    ghost_dot(&rnorm,r,r,A->context->mpicomm);
    rnorm = sqrt(rnorm);
    if(myrank==printRank) 
		{
			printf("|Ax-b|      = %g\n",rnorm);
    	printf("-------------------------------------\n");
   	} 
//    essexamples_print_info(A,0);

    ghost_densemat_destroy(v);
    ghost_densemat_destroy(r);
    ghost_densemat_destroy(x);
    ghost_densemat_destroy(p);
    ghost_densemat_destroy(b);
    ghost_sparsemat_destroy(A);

    ghost_context_destroy(context);
#ifdef AFT	
	AFT_END();
#endif
	
    return NULL;
    
}

int main(int argc, char* argv[])
{
	args arg;
	arg.argc = argc;
	arg.argv = argv;
  essexamples_process_options(argc,argv);
  ghost_init(argc,argv); 
  ghost_task *t;
  ghost_task_create(&t,GHOST_TASK_FILL_ALL,0,&mainTask,&arg,GHOST_TASK_DEFAULT, NULL, 0);

  ghost_task_enqueue(t);
 	ghost_task_wait(t);
  ghost_task_destroy(t);
  ghost_finalize();

  return EXIT_SUCCESS;
}


