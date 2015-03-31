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
                      "graphene<L> or graphene<M>x<N> for an L times L (M times N) graphene sheet,\n"
                      "spinSZ<L> for a spin chain of L spins\n"
                      "anderson<L> for an L^3 model problem for the Anderson localization\n"
                      "matpde<L> for an L^2 eigenproblem from a scalar elliptic partial \n"
                      "          differential equation\n");
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
GRAPHENE,
ANDERSON,
SPINSZ,
MATPDE
} problem_t;

// definitions for MATPDE
void MATPDE_initDimensions(int, int, gidx_t*, lidx_t*);
int MATPDE_rowFunc(gidx_t, lidx_t*, gidx_t*, void*);

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

  int L; // problem size for Graphene (L x L grid) or the Anderson model (L^3 grid)

  if (strcmp(problem,"usage")==0)
  {
    *iflag=0;
    SUBR(create_matrix_usage)();
    return;
  }

  // check if the first arg is of the grapheneN or andersonN form,
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
    mat_type=ANDERSON;
    pos=strlen("anderson");
  }
  else if( str_starts_with(problem,"matpde") )
  {
    mat_type=MATPDE;
    pos=strlen("matpde");
  }

  gidx_t DIM;

  // if not from file and not from bapps, it is one
  // of the matfuncs cases. Check that the string is
  // correctly formatted, otherwise assume FROM_FILE
  // after all.
  if (mat_type!=FROM_FILE && mat_type!=FROM_BAPPS)
  {
    // make sure all remaining characters are
    for (int i=pos;i<strlen(problem);i++)
    {
      if (problem[i]<'0' || problem[i]>'9' && problem[i]!='x')
      {
        mat_type=FROM_FILE;
        break;
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
    PHIST_SOUT(PHIST_INFO,"problem type: Graphene %d x %d\n",L1,L2);
    double gamma=0.2;
    ghost_lidx_t WL[2];
    WL[0] = L1;
    WL[1] = L2;
  
    matfuncs_info_t info;
    // set problem size
    crsGraphene( -2, WL, &DIM, NULL);
    // set disorder
    crsGraphene( -5, NULL, NULL,&gamma);
    // get matrix info
    crsGraphene( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &crsGraphene, iflag), *iflag);
  }
  else if (mat_type==ANDERSON)
  {
    ghost_lidx_t LL[3];
    LL[0]=L,LL[1]=L,LL[2]=L;

    int L1=L;
    int L2=L;
    for (int j=1;j<=2;j++) 
    {
      for (int i=pos; i<strlen(problem);i++)
      {
        if (problem[i]=='x')
        {
          pos=i+1;
          break;
        }
      }
      LL[j]=atoi(problem+pos);
      if (LL[j]<0)
      {
        LL[j]=L;
        break;
      }
    }
  
    PHIST_SOUT(PHIST_INFO,"problem type: Anderson %d x %d x %d\n",LL[0],LL[1],LL[2]);
#ifdef PHIST_HAVE_ESSEX_PHYSICS
        PHIST_SOUT(PHIST_ERROR,"ERROR anderson model problem disabled with ESSEX-Physics\n" 
                               "due to Issue #107\n");
        *iflag=-1;
        return;
#else
    matfuncs_info_t info;
    anderson( -2, LL, &DIM, NULL);
    anderson( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &anderson, iflag), *iflag);
#endif
  }
  else if (mat_type==SPINSZ)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: spinSZ[%d]\n",L);

    ghost_lidx_t conf_spinZ[3] = {L,L/2,0};
    SpinChainSZ( -2, conf_spinZ, &DIM, NULL);

    matfuncs_info_t info;
    SpinChainSZ( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &SpinChainSZ, iflag), *iflag);
  }
  else if (mat_type==MATPDE)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: MATPDE %d x %d\n", L, L);
    gidx_t nrows = -1;
    gidx_t ncols = -1;
    lidx_t row_nnz = -1;
    MATPDE_initDimensions(L, L, &nrows, &row_nnz);
    ncols = nrows;
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &MATPDE_rowFunc, iflag), *iflag);
  }
  else if (mat_type==FROM_BAPPS)
  {
#ifndef PHIST_HAVE_ESSEX_PHYSICS
PHIST_SOUT(PHIST_ERROR,"BAPPS models (essex-physics/bapps) not\n"
                       "available without essex-physics\n");
  *iflag=-99;
#else
    
    int len=strlen(problem+pos);
    // we need to copy the string because bapp_select_model
    // uses strsep, which destroys the string.
    char *matstr=(char*)malloc(len+1);
    strncpy(matstr,problem+pos,len+1);
    PHIST_SOUT(PHIST_DEBUG,"problem '%s', matstr '%s'\n",problem,matstr);
    ghost_sparsemat_fromRowFunc_t matfunc;
    bapp_select_model(matstr,&matfunc);
    free(matstr);
    
    matfuncs_info_t info;
                 
    // query
    matfunc( -2, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFunc)(mat, comm, 
          info.nrows, info.nrows, info.row_nnz, matfunc, iflag), *iflag);  
#endif
  }
  else
  {
    PHIST_SOUT(PHIST_INFO,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(sparseMat_read)(mat,comm,(char*)problem,iflag),*iflag);
  }
  
  return;
}
#else
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* iflag)
{
  *iflag=0;
  if (strcmp(problem,"usage")==0)
  {
    *iflag=0;
    SUBR(create_matrix_usage)();
    PHIST_SOUT(PHIST_INFO,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(sparseMat_read)(mat,comm,(char*)problem,iflag),*iflag);
    return;
  }

}
#endif
