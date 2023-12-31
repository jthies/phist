/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
// auto-detects the file type by looking at the file extension. Either comm or map must be non-NULL
void SUBR(sparseMat_read)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        phist_const_context_ptr ctx, char* filename, int* iflag)
{
  if (phist_filename_isMM(filename))
  {
    if (ctx)
    {
      PHIST_CHK_IERR(SUBR(sparseMat_read_mm_with_context)(A,ctx,filename,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(sparseMat_read_mm)(A,comm,filename,iflag),*iflag);
    }
  }
  else if (phist_filename_isHB(filename))
  {
    if (ctx)
    {
      PHIST_CHK_IERR(SUBR(sparseMat_read_hb_with_context)(A,ctx,filename,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(sparseMat_read_hb)(A,comm,filename,iflag),*iflag);
    }
  }
  else if (phist_filename_isCRS(filename))
  {
    if (ctx)
    {
      PHIST_CHK_IERR(SUBR(sparseMat_read_bin_with_context)(A,ctx,filename,iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(SUBR(sparseMat_read_bin)(A,comm,filename,iflag),*iflag);
    }
  }
  else
  {
      *iflag=-1;
  }
  return;
}

void SUBR(create_matrix_usage)(void)
{
    int iflag;
#ifdef PHIST_HAVE_MPI
    MPI_Comm comm=MPI_COMM_WORLD;
#else
    MPI_Comm comm = 0;
#endif
    PHIST_SOUT(PHIST_INFO,"\nsupported matrix formats:\n\n");
    SUBR(sparseMat_read_mm)(NULL,&comm,NULL,&iflag);
    if (iflag!=PHIST_NOT_IMPLEMENTED)
    {
      PHIST_SOUT(PHIST_INFO," * Matrix Market ascii format (.mm|.mtx)\n");
    }
    SUBR(sparseMat_read_hb)(NULL,&comm,NULL,&iflag);
    if (iflag!=PHIST_NOT_IMPLEMENTED)
    {
#ifdef IS_COMPLEX
      PHIST_SOUT(PHIST_INFO," * Harwell Boeing ascii format (.cua)\n");
#else
      PHIST_SOUT(PHIST_INFO," * Harwell Boeing ascii format (.rua)\n");
#endif
    }
    SUBR(sparseMat_read_bin)(NULL,&comm,NULL,&iflag);
    if (iflag!=PHIST_NOT_IMPLEMENTED)
    {
      PHIST_SOUT(PHIST_INFO," * GHOST binary CRS format (.bin|.crs)\n");
    }

#if defined(IS_DOUBLE) && !defined(IS_COMPLEX)
PHIST_SOUT(PHIST_INFO,"\n\nInstead of a matrix file you can also specify a string describing\n"
                      "one of our scalable test problems, e.g. \n"
                      "graphene<L> or graphene<M>x<N> for an L times L (M times N) graphene sheet\n"
                      "spin<L> for a spin chain of L spins\n"
                      "spinSZ<L> for a spin chain of L spins with SZ-symmetry\n"
                      "matpde<L> for an L^2 eigenproblem from a scalar elliptic partial differential equation\n"
                      "TriToeplitz<L> for an 2^L tridiagonal Toeplitz matrix (e.g. 1D Poisson, spd, diagonal dominant)\n"
                      "Brussolator<L> for an L^1 eigenproblem from a Brussolator wave model in chemical reaction (MVMBWM)\n"
                      "BENCH3D-<L>-<WHICH> for an L^3 scalar problem. \n"
                      "                    the strings have the form <A,B,C><0-9>, for instance\n"
                      "                    B1 is an inhomogenous PDE problem similar to matpde but in 3D, C1 is the Anderson model problem\n"
                      "                    and the 'A' problems are convection diffusion problems taken from Gordon&Gordon, ParCo'10.\n"
                      "                    The '0' cases are symmetric matrices, A0=Laplace, B0=variable coefficient diffusion.\n"
# ifdef PHIST_HAVE_SCAMAC
                      "SCAMAC-<string> adds a variety of benchmarks from the ESSEX SCAlable MAtrix Collection, for example\n"
                      "                    'SCAMAC-hubbard,n_sites=10' etc.\n"
                      "                Use the driver 'scamact' to find out more about the options.\n"
                      "                The stand-alone repo of the scamac can be found at https://bitbucket.org/essex/MatrixCollection\n"
# endif
                      );
#endif

}

// variant for double real that allows creating selected matrices from functions
#if defined(IS_DOUBLE) && !defined(IS_COMPLEX)

// ghost/physics
#ifdef PHIST_KERNEL_LIB_GHOST
#include <ghost.h>
#endif
#include "matfuncs.h"
#ifdef PHIST_HAVE_SCAMAC
#include "scamac.h"
#endif

typedef enum {
FROM_FILE,
FROM_SCAMAC,
FROM_BENCH3D,
GRAPHENE,
SPIN,
SPINSZ,
MATPDE,
TRITOEPLITZ,
BRUSSOLATOR
} problem_t;

// definitions for MATPDE
void MATPDE_initDimensions(int, int, phist_gidx*, phist_lidx*);
int MATPDE_rowFunc(phist_gidx, phist_lidx*, phist_gidx*, void*, void*);
// definitions for TriToeplitz
void TriToeplitz_initDimensions(int, phist_gidx*, phist_lidx*);
int TriToeplitz_rowFunc(phist_gidx, phist_lidx*, phist_gidx*, void*, void*);
// definitions for MATPDE3D
void MATPDE3D_initDimensions(int, int, int, phist_gidx*, phist_lidx*);
void MATPDE3D_selectProblem(int which, int* iflag);
int MATPDE3D_rowFunc(phist_gidx, phist_lidx*, phist_gidx*, void*, void*);
int MATPDE3D_solFunc(phist_gidx, phist_lidx, void*,void*);
int MATPDE3D_rhsFunc(phist_gidx, phist_lidx, void*,void*);
// definitions for Brussolator
void Brussolator_initDimensions(int, phist_gidx*, phist_lidx*);
int Brussolator_rowFunc(phist_gidx, phist_lidx*, phist_gidx*, void*, void*);

int str_starts_with(const char *s1, const char *s2)
{
  int minLen = strlen(s2);
  if( strlen(s1) < minLen )
    return 0;

  if( strncmp(s1, s2, minLen) == 0 )
    return 1;

  return 0;
}

#ifdef PHIST_HAVE_SCAMAC

//////////////////////////////////////////////////////////////////////////
// helper struct and functions for the scamac (ESSEX matrix collection) //
// copied from ghost/minimal_scamac example by Andreas Alvermann        //
//////////////////////////////////////////////////////////////////////////

struct phist_scamac_work_st {
 ScamacGenerator * gen;
 ScamacWorkspace * ws;
};

static int phist_scamac_func(ghost_gidx row, ghost_lidx *rowlen, ghost_gidx *col, void *val, void *arg)
{
    struct phist_scamac_work_st * my_work = arg;
    ScamacIdx nzr;
    ScamacErrorCode err=scamac_generate_row(my_work->gen, my_work->ws, (ScamacIdx) row, SCAMAC_DEFAULT, &nzr, (ScamacIdx *) col, (double *) val); 
    *rowlen=(ghost_lidx)nzr;
    if (err) { return 1; }
    return 0;
}

static int phist_scamac_funcinit(void *arg, void **work) 
{
    int iflag=0;
    PHIST_ICHK_IERR(iflag=(sizeof(ScamacIdx)==sizeof(ghost_gidx))?0:-1,iflag);
  ScamacGenerator * gen = (ScamacGenerator *) arg;
  if (*work) {// free
    struct phist_scamac_work_st * my_work = *work;
    scamac_workspace_free( my_work->ws);
    free(my_work);
    *work=NULL;
  } else {//alloc
    struct phist_scamac_work_st * my_work = malloc(sizeof *my_work);
    my_work->gen = (ScamacGenerator * ) arg;
    scamac_workspace_alloc(my_work->gen, &(my_work->ws));
    *work = my_work;
  } 
  return iflag;
}

#endif

// TODO: we should do this in C++ with a map where you can easily add new function pointers... (melven)
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, phist_const_comm_ptr comm,
        const char* problem, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  problem_t mat_type=FROM_FILE; // default: assume that 'which' is a file name
  int outlev = (*iflag & PHIST_SPARSEMAT_QUIET) ? PHIST_DEBUG : PHIST_VERBOSE;

  int L; // problem size for Graphene (L x L grid) or the Anderson model (L^3 grid)

  if (strcmp(problem,"usage")==0)
  {
    *iflag=0;
    SUBR(create_matrix_usage)();
    return;
  }

  // check if the first arg is a supported matrix string
  // otherwise try to read a file
  mat_type=FROM_FILE;
  int pos=0;
  if( str_starts_with(problem,"SCAMAC-") )
  {
#ifdef PHIST_HAVE_SCAMAC
    mat_type=FROM_SCAMAC;
    pos=strlen("SCAMAC-");
#else
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
#endif
  }
  else if( str_starts_with(problem,"spinSZ") )
  {
    mat_type=SPINSZ;
    pos=strlen("spinSZ");
  }
  else if( str_starts_with(problem,"spin") )
  {
    mat_type=SPIN;
    pos=strlen("spin");
  }
  else if( str_starts_with(problem,"graphene") )
  {
    mat_type=GRAPHENE;
    pos=strlen("graphene");
  }
  else if( str_starts_with(problem,"anderson") )
  {
    PHIST_SOUT(PHIST_WARNING,"Problem 'anderson<L>' is now generated by 'BENCH3D-<L>-C0', \n"
                             "assuming the string you passed to create_matrix is a file name.\n"
                             "(file %s, line %d)\n",__FILE__,__LINE__);
  }
  else if( str_starts_with(problem,"BENCH3D-") )
  {
    mat_type=FROM_BENCH3D;
    pos=strlen("BENCH3D-");
  }
  else if( str_starts_with(problem,"matpde") )
  {
    mat_type=MATPDE;
    pos=strlen("matpde");
  }
  else if( str_starts_with(problem,"TriToeplitz") )
  {
    mat_type=TRITOEPLITZ;
    pos=strlen("TriToeplitz");
  }
  else if( str_starts_with(problem,"Brussolator") )
  {
    mat_type=BRUSSOLATOR;
    pos=strlen("Brussolator");
  }


  phist_gidx DIM;

  // if not from file and not from bapps, it is one
  // of the matfuncs cases. Check that the string is
  // correctly formatted, otherwise assume FROM_FILE
  // after all.
  if (mat_type!=FROM_FILE && mat_type!=FROM_SCAMAC)
  {
    // make sure all remaining characters are
    if (mat_type!=FROM_BENCH3D)
    {
      for (int i=pos;i<strlen(problem);i++)
      {
        if ( (problem[i]<'0' || problem[i]>'9') && problem[i]!='x')
        {
          mat_type=FROM_FILE;
          break;
        }
      }
    }

    const char* strL=problem+pos;
    L=atoi(strL);
    if (L<0) mat_type=FROM_FILE;
  }

  // check for matfuncs cases
  if (mat_type==GRAPHENE)
  {
    int L1=L;
    int L2=L;
    for (int i=pos; i<strlen(problem);i++)
    {
      if (problem[i]=='x')
      {
        pos=i+1;
        break;
      }
    }
    L2=atoi(problem+pos);
    if (L2<0)
    {
      L2=L1;
    }
    PHIST_SOUT(outlev,"problem type: Graphene %d x %d\n",L1,L2);
    double gamma=0.2;
    ghost_lidx WL[2];
    WL[0] = L1;
    WL[1] = L2;
  
    matfuncs_info_t info;
    // set problem size
    crsGraphene( -2, WL, &DIM, NULL, NULL);
    // set disorder
    crsGraphene( -5, NULL, NULL,&gamma, NULL);
    // get matrix info
    crsGraphene( -1, NULL, NULL, &info, NULL);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (phist_gidx)info.nrows, (phist_gidx)info.ncols, (phist_lidx)info.row_nnz,
        &crsGraphene, NULL, iflag), *iflag);
  }
  else if (mat_type==SPIN)
  {
  //TODO!!!
    PHIST_SOUT(outlev,"problem type: spin[%d]\n",L);

    ghost_lidx conf_spin[2] = {L,0};
    SpinChain( -2, conf_spin, &DIM, NULL, NULL);

    matfuncs_info_t info;
    SpinChain( -1, NULL, NULL, &info, NULL);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (phist_gidx)info.nrows, (phist_gidx)info.ncols, (phist_lidx)info.row_nnz,
        &SpinChain, NULL, iflag), *iflag);
  }
  else if (mat_type==SPINSZ)
  {
    PHIST_SOUT(outlev,"problem type: spinSZ[%d]\n",L);

    ghost_lidx conf_spinZ[3] = {L,L/2,0};
    SpinChainSZ( -2, conf_spinZ, &DIM, NULL, NULL);

    matfuncs_info_t info;
    SpinChainSZ( -1, NULL, NULL, &info, NULL);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (phist_gidx)info.nrows, (phist_gidx)info.ncols, (phist_lidx)info.row_nnz,
        &SpinChainSZ, NULL, iflag), *iflag);
    // clean up internal memory allocated by the function
    SpinChainSZ( -3, conf_spinZ, &DIM, NULL, NULL);
  }
  else if (mat_type==MATPDE)
  {
    PHIST_SOUT(outlev,"problem type: MATPDE %d x %d\n", L, L);

    phist_gidx nrows = -1;
    phist_lidx row_nnz = -1;
    MATPDE_initDimensions(L, L, &nrows, &row_nnz);
    phist_gidx ncols = nrows;
    if( *iflag & PHIST_SPARSEMAT_PERM_GLOBAL )
    {
      PHIST_SOUT(outlev,"Disabling PHIST_SPARSEMAT_PERM_GLOBAL; MATPDE features a predefined partitioning!\n");
      *iflag &= ~PHIST_SPARSEMAT_PERM_GLOBAL;
    }
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &MATPDE_rowFunc, NULL, iflag), *iflag);
  }
  else if(mat_type==TRITOEPLITZ)
  {
    PHIST_SOUT(outlev,"problem type: TriToeplitz 2^%d\n", L);
    phist_gidx nrows = -1;
    phist_lidx row_nnz = -1;
    TriToeplitz_initDimensions(L, &nrows, &row_nnz);
    phist_gidx ncols = nrows;
    if( *iflag & PHIST_SPARSEMAT_PERM_GLOBAL )
    {
      PHIST_SOUT(outlev,"Disabling PHIST_SPARSEMAT_PERM_GLOBAL; TriToeplitz is already tridiagonal!\n");
      *iflag &= ~PHIST_SPARSEMAT_PERM_GLOBAL;
    }
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &TriToeplitz_rowFunc, NULL, iflag), *iflag);
  }
  else if(mat_type==BRUSSOLATOR)
  {
    PHIST_SOUT(outlev,"problem type: Brussolator wave model %d^1\n", L);
    phist_gidx nrows = -1;
    phist_lidx row_nnz = -1;
    Brussolator_initDimensions(L, &nrows, &row_nnz);
    phist_gidx ncols = nrows;
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &Brussolator_rowFunc, NULL, iflag), *iflag);
  }
  else if (mat_type==FROM_BENCH3D)
  {
      for (int i=pos; i<strlen(problem);i++)
      {
        if (problem[i]=='-')
        {
          pos=i+1;
          break;
        }
      }

    long int which=strtol(problem+pos, NULL, 16);

    PHIST_SOUT(outlev,"problem type: BENCH3D %x %dx%dx%d\n", (uint32_t)which,L, L, L);

    phist_gidx nrows = -1;
    phist_lidx row_nnz = -1;
    MATPDE3D_initDimensions(L, L, L, &nrows, &row_nnz);
    phist_gidx ncols = nrows;
    if( *iflag & PHIST_SPARSEMAT_PERM_GLOBAL )
    {
      PHIST_SOUT(outlev,"Disabling PHIST_SPARSEMAT_PERM_GLOBAL; MATPDE3D features a predefined partitioning!\n");
      *iflag &= ~PHIST_SPARSEMAT_PERM_GLOBAL;
    }
    int iflag_tmp=*iflag;
    PHIST_CHK_IERR(MATPDE3D_selectProblem(which,iflag),*iflag);
    *iflag=iflag_tmp;
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &MATPDE3D_rowFunc, NULL, iflag), *iflag);
  }
  else if (mat_type==FROM_SCAMAC)
  {
#ifdef PHIST_HAVE_SCAMAC
    int len=strlen(problem+pos);
    // we need to copy the string because bapp_select_model
    // uses strsep, which destroys the string.
    char *matstr=(char*)malloc(len+1);
    strncpy(matstr,problem+pos,len+1);
    PHIST_SOUT(outlev,"problem '%s', matstr '%s'\n",problem,matstr);
    ScamacGenerator* gen;
    char* errdesc=NULL;
    ScamacErrorCode scamac_err=scamac_parse_argstr(matstr, &gen, &errdesc);
    free(matstr);
    
    if (scamac_err!=SCAMAC_EOK)
    {
      PHIST_SOUT(PHIST_ERROR,"%s\n",errdesc);
      PHIST_CHK_IERR(*iflag=(int)scamac_err, *iflag);
    }

    PHIST_CHK_IERR(*iflag=scamac_err=scamac_generator_finalize(gen),*iflag);
    
    // query matrix size and max num nonzereos per row
    ghost_gidx gnrows = scamac_generator_query_nrow(gen);
    int max_nnz = scamac_generator_query_maxnzrow(gen);


    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncWithConstructor)(mat, comm, 
          gnrows, gnrows, max_nnz, phist_scamac_func, phist_scamac_funcinit, gen, iflag), *iflag);  
    scamac_generator_destroy(gen);  
#else
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
#endif
  }
  else
  {
    PHIST_SOUT(outlev,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(sparseMat_read)(mat,comm,NULL,(char*)problem,iflag),*iflag);
  }
  
  return;
}
#else
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, phist_const_comm_ptr comm,
        const char* problem, int* iflag)
{
  if (strcmp(problem,"usage")==0)
  {
    SUBR(create_matrix_usage)();
    return;
  }

    PHIST_SOUT(PHIST_INFO,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(sparseMat_read)(mat,comm,NULL,(char*)problem,iflag),*iflag);

}
#endif

void SUBR(create_matrix_with_context)(TYPE(sparseMat_ptr)* mat, phist_const_context_ptr ctx,
        const char* problem, int* iflag)
{
  // try interpreting the string as filename
  SUBR(sparseMat_read)(mat,NULL,ctx,(char*)problem,iflag);
  // if that fails, return NOT_IMPLEMENTED (input from a function can be added later if we need it).
  // Note that we do have a kernel function sparseMat_create_fromRowFuncAndMap, but interpreting the
  // string to select and initialize the rowFunc is currently done in the create_matrix function.
  if (*iflag!=0) PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
}

//! For testing linear solvers, generates an 'exact solution' sol and right-hand side rhs
//! for some matrix creted by create_matrix. For most test cases, this will be some random
//! sol vector and the rhs is computed by rhs=A*sol, but for the cases stemming from PDEs
//! (BENCH3D-A*,B*) we presribe an analytical solution and generate the correct F for it.
//! The mvecs sol and rhs must be created beforehand and may have an arbitrary number of 
//! columns.
void SUBR(create_sol_and_rhs)(const char* problem, TYPE(const_sparseMat_ptr) A,
                        TYPE(mvec_ptr) sol, TYPE(mvec_ptr) rhs, int* iflag)
{

  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  int sol_done=0, rhs_done=0;
  // only the MATPDE3D Benchmark problems have
  // an implementation for providing an analytical
  // solution and exact RHS vector. We skip the rhs
  // though because the analytical rhs takes discretization
  // errors into account when used for residual calculations.
  // So we use the analytical solution x and then compute b=A*x.
  // For the other
  // cases (and matrices read from a file), we use
  // a random vector X and compute B=A*X.
#if defined(IS_DOUBLE) && !defined(IS_COMPLEX)
  if (str_starts_with(problem,"BENCH3D"))
  {
    // MATPDE3D provides these for us
    SUBR(mvec_put_func)(sol,&MATPDE3D_solFunc,NULL,iflag);
    if (*iflag==0)
    {
      sol_done=1;
    }
//    if (sol_done)
//    {
//      SUBR(mvec_put_func)(rhs,&MATPDE3D_rhsFunc,NULL,iflag);
//    if (*iflag==0) rhs_done=1;
//    }
  }
#endif
  if (!sol_done)
  {
    PHIST_CHK_IERR(SUBR(mvec_random)(sol,iflag),*iflag);
  }
  if (!rhs_done)
  {
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(ONE,A,sol,ZERO,rhs,iflag),*iflag);
  }
}
