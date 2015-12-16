// auto-detects the file type by looking at the file extension
void SUBR(sparseMat_read)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm,
        char* filename, int* iflag)
{
  char *c=filename+strlen(filename)-1;
  int isMM=0;
  int isCRS=0;
  int isHB=0;
  
  while (c!=filename)
  {
    if (*c=='.') break;
    c--;
  }
  c++;
  isMM=(!strcmp(c,"mm"))||(!strcmp(c,"mtx"));
  isCRS=!strcmp(c,"crs");
  isCRS=isCRS||!strcmp(c,"bin");
  isHB=!strcmp(c,"rua")||!strcmp(c,"cua");
  if (isMM)
  {
    PHIST_CHK_IERR(SUBR(sparseMat_read_mm)(A,comm,filename,iflag),*iflag);
  }
  else if (isHB)
  {
    PHIST_CHK_IERR(SUBR(sparseMat_read_hb)(A,comm,filename,iflag),*iflag);
  }
  else if (isCRS)
  {
    PHIST_CHK_IERR(SUBR(sparseMat_read_bin)(A,comm,filename,iflag),*iflag);
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
    MPI_Comm comm=MPI_COMM_WORLD;
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
                      "spinSZ<L> for a spin chain of L spins\n"
                      "matpde<L> for an L^2 eigenproblem from a scalar elliptic partial differential equation\n"
                      "TriToeplitz<L> for an 2^L tridiagonal Toeplitz matrix (e.g. 1D Poisson, spd, diagonal dominant)\n"
                      "Brussolator<L> for an L^1 eigenproblem from a Brussolator wave model in chemical reaction (MVMBWM)\n"
                      "BENCH3D-<L>-<WHICH> for an L^3 scalar problem. \n"
                      "                    the strings have the form <A,B,C><0-9>, for instance\n"
                      "                    B1 is an inhomogenous PDE problem similar to matpde but in 3D, C1 is the Anderson model problem\n"
                      "                    and the 'A' problems are convection diffusion problems taken from Gordon&Gordon, ParCo'10.\n"
                      "                    The '0' cases are symmetric matrices, A0=Laplace, B0=variable coefficient diffusion.\n"
                      "BAPP-<string> requires the ESSEX-Physics library and adds a variety of benchmarks, for example\n"
                      "                    'BAPP-hubbard,l=10,U=1' etc.\n"
                      );
#endif

}

#if defined(IS_DOUBLE) && !defined(IS_COMPLEX)

// ghost/physics
#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
#endif
#ifdef PHIST_HAVE_ESSEX_PHYSICS
#include "essex-physics/matfuncs.h"
#include "essex-physics/bapp.h"
#else
#include "matfuncs.h"
#endif
typedef enum {
FROM_FILE,
FROM_BAPPS,
FROM_BENCH3D,
GRAPHENE,
SPINSZ,
MATPDE,
TRITOEPLITZ,
BRUSSOLATOR
} problem_t;

// definitions for MATPDE
void MATPDE_initDimensions(int, int, gidx_t*, lidx_t*);
int MATPDE_rowFunc(gidx_t, lidx_t*, gidx_t*, void*, void*);
// definitions for TriToeplitz
void TriToeplitz_initDimensions(int, gidx_t*, lidx_t*);
int TriToeplitz_rowFunc(gidx_t, lidx_t*, gidx_t*, void*, void*);
// definitions for MATPDE3D
void MATPDE3D_initDimensions(int, int, int, gidx_t*, lidx_t*);
void MATPDE3D_selectProblem(int which, int* iflag);
int MATPDE3D_rowFunc(gidx_t, lidx_t*, gidx_t*, void*, void*);
int MATPDE3D_solFunc(gidx_t, lidx_t, void*,void*);
int MATPDE3D_rhsFunc(gidx_t, lidx_t, void*,void*);
// definitions for Brussolator
void Brussolator_initDimensions(int, gidx_t*, lidx_t*);
int Brussolator_rowFunc(gidx_t, lidx_t*, gidx_t*, void*, void*);

int str_starts_with(const char *s1, const char *s2)
{
  int minLen = strlen(s2);
  if( strlen(s1) < minLen )
    return 0;

  if( strncmp(s1, s2, minLen) == 0 )
    return 1;

  return 0;
}

// TODO: we should do this in C++ with a map where you can easily add new function pointers... (melven)
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  problem_t mat_type=FROM_FILE; // default: assume that 'which' is a file name
  int outlev = *iflag & PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;
  outlev = PHIST_INFO;

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
  if( str_starts_with(problem,"BAPP-") )
  {
    mat_type=FROM_BAPPS;
    pos=strlen("BAPP-");
  }
  else if( str_starts_with(problem,"spinSZ") )
  {
    mat_type=SPINSZ;
    pos=strlen("spinSZ");
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


  gidx_t DIM;

  // if not from file and not from bapps, it is one
  // of the matfuncs cases. Check that the string is
  // correctly formatted, otherwise assume FROM_FILE
  // after all.
  if (mat_type!=FROM_FILE && mat_type!=FROM_BAPPS)
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
    ghost_lidx_t WL[2];
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
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &crsGraphene, NULL, iflag), *iflag);
  }
  else if (mat_type==SPINSZ)
  {
    PHIST_SOUT(outlev,"problem type: spinSZ[%d]\n",L);

    ghost_lidx_t conf_spinZ[3] = {L,L/2,0};
    SpinChainSZ( -2, conf_spinZ, &DIM, NULL, NULL);

    matfuncs_info_t info;
    SpinChainSZ( -1, NULL, NULL, &info, NULL);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &SpinChainSZ, NULL, iflag), *iflag);
  }
  else if (mat_type==MATPDE)
  {
    PHIST_SOUT(outlev,"problem type: MATPDE %d x %d\n", L, L);

    gidx_t nrows = -1;
    lidx_t row_nnz = -1;
    MATPDE_initDimensions(L, L, &nrows, &row_nnz);
    gidx_t ncols = nrows;
    if( *iflag & PHIST_SPARSEMAT_REPARTITION )
    {
      PHIST_SOUT(outlev,"Disabling PHIST_SPARSEMAT_REPARTITION; MATPDE features a predefined partitioning!\n");
      *iflag &= ~PHIST_SPARSEMAT_REPARTITION;
    }
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &MATPDE_rowFunc, NULL, iflag), *iflag);
  }
  else if(mat_type==TRITOEPLITZ)
  {
    PHIST_SOUT(outlev,"problem type: TriToeplitz 2^%d\n", L);
    gidx_t nrows = -1;
    lidx_t row_nnz = -1;
    TriToeplitz_initDimensions(L, &nrows, &row_nnz);
    gidx_t ncols = nrows;
    if( *iflag & PHIST_SPARSEMAT_REPARTITION )
    {
      PHIST_SOUT(outlev,"Disabling PHIST_SPARSEMAT_REPARTITION; TriToeplitz is already tridiagonal!\n");
      *iflag &= ~PHIST_SPARSEMAT_REPARTITION;
    }
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &TriToeplitz_rowFunc, NULL, iflag), *iflag);
  }
  else if(mat_type==BRUSSOLATOR)
  {
    PHIST_SOUT(outlev,"problem type: Brussolator wave model %d^1\n", L);
    gidx_t nrows = -1;
    lidx_t row_nnz = -1;
    Brussolator_initDimensions(L, &nrows, &row_nnz);
    gidx_t ncols = nrows;
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

    gidx_t nrows = -1;
    lidx_t row_nnz = -1;
    MATPDE3D_initDimensions(L, L, L, &nrows, &row_nnz);
    gidx_t ncols = nrows;
    if( *iflag & PHIST_SPARSEMAT_REPARTITION )
    {
      PHIST_SOUT(outlev,"Disabling PHIST_SPARSEMAT_REPARTITION; MATPDE3D features a predefined partitioning!\n");
      *iflag &= ~PHIST_SPARSEMAT_REPARTITION;
    }
    int iflag_tmp=*iflag;
    PHIST_CHK_IERR(MATPDE3D_selectProblem(which,iflag),*iflag);
    *iflag=iflag_tmp;
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &MATPDE3D_rowFunc, NULL, iflag), *iflag);
  }
  else if (mat_type==FROM_BAPPS)
  {
#ifndef PHIST_HAVE_ESSEX_PHYSICS
PHIST_SOUT(PHIST_ERROR,"BAPPS models (essex-physics/bapps) not\n"
                       "available without essex-physics\n");
  *iflag=PHIST_NOT_IMPLEMENTED;
#else
    
    int len=strlen(problem+pos);
    // we need to copy the string because bapp_select_model
    // uses strsep, which destroys the string.
    char *matstr=(char*)malloc(len+1);
    strncpy(matstr,problem+pos,len+1);
    PHIST_SOUT(outlev,"problem '%s', matstr '%s'\n",problem,matstr);
    ghost_sparsemat_fromRowFunc_t matfunc;
    bapp_select_model(matstr,&matfunc);
    free(matstr);
    
    matfuncs_info_t info;
                 
    // query
    matfunc( -2, NULL, NULL, &info, NULL);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          info.nrows, info.nrows, info.row_nnz, matfunc, NULL, iflag), *iflag);  
#endif
  }
  else
  {
    PHIST_SOUT(outlev,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(sparseMat_read)(mat,comm,(char*)problem,iflag),*iflag);
  }
  
  return;
}
#else
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* iflag)
{
  if (strcmp(problem,"usage")==0)
  {
    SUBR(create_matrix_usage)();
    PHIST_SOUT(PHIST_INFO,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(sparseMat_read)(mat,comm,(char*)problem,iflag),*iflag);
    return;
  }
  *iflag=0;

}
#endif

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
  int done=0;
  // only the MATPDE3D Benchmark problems have
  // an implementation for providing an analytical
  // solution and exact RHS vector, for the other
  // cases (and matrices read from a file), we use
  // a random vector X and compute B=A*X.
#if defined(IS_DOUBLE) && !defined(IS_COMPLEX)
  if (str_starts_with(problem,"BENCH3D"))
  {
    // MATPDE3D provides these for us
    SUBR(mvec_put_func)(sol,&MATPDE3D_solFunc,NULL,iflag);
    if (*iflag==0) done=1;
    SUBR(mvec_put_func)(rhs,&MATPDE3D_rhsFunc,NULL,iflag);
    if (*iflag==0) done&=1;
  }
#endif
  if (!done)
  {
    // if not, or not a BENCH3D case, set sol=random and rhs=A*sol
    PHIST_CHK_IERR(SUBR(mvec_random)(sol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(ONE,A,sol,ZERO,rhs,iflag),*iflag);
  }
}
