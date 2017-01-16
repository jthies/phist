
#include <essexamples.h>
#include <ghost.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <numeric>
#define DP

#ifdef DP
GHOST_REGISTER_DT_D(vecdt) // vectors have double values
GHOST_REGISTER_DT_D(matdt) // matrix has double values
#define SQRT(a) sqrt(a)
#define FABS(a) fabs(a)
#define IMTQL1(a,b,c,d) imtql1_(a,b,c,d)
#else
GHOST_REGISTER_DT_S(vecdt) // vectors have float values
GHOST_REGISTER_DT_S(matdt) // matrix has float values
#define SQRT(a) sqrtf(a)
#define FABS(a) fabsf(a)
#define IMTQL1(a,b,c,d) imtql1f_(a,b,c,d)
#endif

#include "lanczos.h"

#include <checkpoint.hpp>
#include <cpOptions.h>
#include <ghost.h>
#include <ghost/types.h>
extern "C"{
#include "essexamples.h"
}
#ifdef AFT
	#include <aft.h>
	#include <aft_macros.h>
#endif

typedef struct args{
	int argc;
	char **argv;
}args;

static int converged(matdt_t evmin)
{
    static matdt_t oldevmin = -1e9;

    int converged = FABS(evmin-oldevmin) < 1e-9;
    oldevmin = evmin;

    return converged;
}

static void switchPtr(ghost_densemat* vnew, ghost_densemat* vold){
  static ghost_densemat* tmp;
  tmp = vnew;
  vnew = vold;
  vold = tmp;
}

static void lanczosStep(ghost_context *context, ghost_sparsemat *mat, ghost_densemat *vnew, ghost_densemat *vold,
        matdt_t *alpha, matdt_t *beta, ghost_spmv_opts spmvtraits)
{
    matdt_t minusbeta = -*beta;
    ghost_scale(vnew,&minusbeta);
    ghost_spmv(vnew, mat, vold, spmvtraits);
    ghost_dot(alpha,vnew,vold);
    matdt_t minusalpha = -*alpha;
    ghost_axpy(vnew,vold,&minusalpha);
    ghost_dot(beta,vnew,vnew);
    *beta=SQRT(*beta);
    matdt_t recbeta = (matdt_t)1./(*beta);
    ghost_scale(vnew,&recbeta);
}

static void *mainTask(void *varg)
{
	args * arg = (args * ) varg;
	int argc = arg->argc;
	char **argv = arg->argv;



    ghost_spmv_opts spmvtraits = GHOST_SPMV_OPTS_INITIALIZER;
    matdt_t alpha=0., beta=0.;
    int ferr, n, iteration= 0, nIter = 500;
    int cp_freq = 999999;
	char * cpPath = new char[256];
	matdt_t zero = 0.;
	matdt_t one = 1.0;
    int myrank;
    const int printrank = 0;


	MPI_Comm FtComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &FtComm);
#ifdef AFT
  AFT_BEGIN(FtComm, &myrank, argv);
#endif 

//	KSM_MAIN_BEGIN(FT_comm, &myrank, argv);

	essexamples_get_iterations(&nIter);
	essexamples_get_cp_freq(&cp_freq);
	essexamples_get_cp_folder(&cpPath);

    ghost_context *context;
    ghost_sparsemat *mat;
    ghost_densemat *vold;
    ghost_densemat *vnew;
    ghost_densemat *tmp;
    
    ghost_sparsemat_traits mtraits = GHOST_SPARSEMAT_TRAITS_INITIALIZER;
    ghost_densemat_traits vtraits = GHOST_DENSEMAT_TRAITS_INITIALIZER;

    essexamples_create_context_and_matrix_ft(&context,&mat,&mtraits, FtComm);
    essexamples_set_spmv_flags(&spmvtraits.flags);
    
    spmvtraits.flags = spmvtraits.flags | GHOST_SPMV_AXPY;

    ghost_rank(&myrank, context->mpicomm);

    ghost_densemat_create(&vnew,ghost_context_max_map(context),vtraits);
    ghost_densemat_create(&vold,ghost_context_max_map(context),vtraits);

    ghost_densemat_init_val(vnew,&zero); // vnew = 0
    ghost_densemat_init_rand(vold); // vold = random
//    ghost_densemat_init_val(vold, &one); 
//    ghost_normalize(vold); // normalize vold    ghost_densemat_permute(vold,GHOST_PERMUTATION_ORIG2PERM);
//    ghost_densemat_permute(vold,GHOST_PERMUTATION_ORIG2PERM);

    matdt_t *alphas  = (matdt_t *)malloc(sizeof(matdt_t)*nIter);
    matdt_t *betas   = (matdt_t *)malloc(sizeof(matdt_t)*nIter);
    matdt_t *falphas = (matdt_t *)malloc(sizeof(matdt_t)*nIter);
    matdt_t *fbetas  = (matdt_t *)malloc(sizeof(matdt_t)*nIter);

    if (!alphas || !betas || !falphas || !fbetas) {
        fprintf(stderr,"Error in malloc!\n");
        exit(EXIT_FAILURE);
    }

    betas[0] = beta;

	int tmp_len = 5; 
	double * a_tmp = new double[tmp_len];
	double * b_tmp = new double[tmp_len];
	for (int i = 0; i < tmp_len ; ++i){
		a_tmp[i] = 0.0;
		b_tmp[i] = 100.0;
	}

	Checkpoint myCP("a", MPI_COMM_WORLD);
	myCP.add("iteration", &iteration);
	myCP.add("alpha", &alpha);
	myCP.add("beta", &beta );
	myCP.add("alphas", alphas, nIter);
	myCP.add("betas", betas, nIter);
	myCP.add("vold", vold);
	myCP.add("vnew", vnew);
	myCP.commit();

	int take_cp = true; 
	int firstrun = true;
	
	if(myCP.needRestart()){
		printf("RESTART ----> \n");
		myCP.read();
		iteration  += 1;
    switchPtr(vnew, vold);  // this pointer-switching is necessary. as vnew and vold are stored at the end of iteration and for each iteration, the pointer has been switched before as well.
	}
    
	for(iteration = 0; iteration < nIter && !converged(falphas[0]); iteration++) 
  {
		n = iteration +1;
    lanczosStep(context,mat,vnew,vold,&alpha,&beta,spmvtraits);
    switchPtr(vnew, vold);

    alphas[iteration] = alpha;
    betas[iteration+1] = beta;
    memcpy(falphas,alphas,nIter*sizeof(matdt_t)); // alphas and betas will be destroyed in imtql
    memcpy(fbetas,betas,nIter*sizeof(matdt_t));

//    IMTQL1(&n,falphas,fbetas,&ferr);

    if(ferr != 0) printf("Error: the %d. eigenvalue could not be determined\n",ferr);
    if (myrank == printrank) {
      printf("%d: min/max eigenvalue: %f/%f\n", iteration, falphas[0],falphas[n-1]);
    }

		if(iteration % cp_freq == 0){
			myCP.update();
		  myCP.write();
		}
  }
    if (myrank == printrank) {
        printf("%s, iterations: %d\n",converged(falphas[0])?" (converged!)":" (max. iterations reached!)",iteration);
    }
//    essexamples_print_info(mat,printrank);

    ghost_densemat_destroy(vold);
    ghost_densemat_destroy(vnew);
    ghost_context_destroy(context);
    ghost_sparsemat_destroy(mat);

#ifdef AFT
	AFT_END()
#endif


//	KSM_MAIN_END();
    return NULL;

}

int main(int argc, char* argv[])
{
	args arg;
	arg.argc = argc;
	arg.argv = argv;

    essexamples_process_options(argc,argv);
    ghost_init(argc,argv); // has to be the first call

    ghost_task *t;
    ghost_task_create(&t,GHOST_TASK_FILL_ALL,0,&mainTask,&arg,GHOST_TASK_DEFAULT, NULL, 0);
    ghost_task_enqueue(t);
    ghost_task_wait(t);
    ghost_task_destroy(t);

    ghost_finalize();
    
    return EXIT_SUCCESS;
}

