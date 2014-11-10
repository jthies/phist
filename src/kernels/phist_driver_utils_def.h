// auto-detects the file type by looking at the file extension
void SUBR(crsMat_read)(TYPE(crsMat_ptr)* A, const_comm_ptr_t comm,
        char* filename, int* ierr)
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
    PHIST_CHK_IERR(SUBR(crsMat_read_mm)(A,comm,filename,ierr),*ierr);
  }
  else if (isHB)
  {
    PHIST_CHK_IERR(SUBR(crsMat_read_hb)(A,comm,filename,ierr),*ierr);
  }
  else if (isCRS)
  {
    PHIST_CHK_IERR(SUBR(crsMat_read_bin)(A,comm,filename,ierr),*ierr);
  }
  else
  {
      *ierr=-1;
  }
  return;
}

#if defined(IS_DOUBLE) && !defined(IS_COMPLEX)

// ghost/physics
#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
#endif
#ifdef PHIST_HAVE_ESSEX_PHYSICS
#include "essex-physics/matfuncs.h"
#else
#include "matfuncs.h"
#endif
typedef enum {
FROM_FILE=0,
GRAPHENE=1,
ANDERSON=2,
SPINSZ=3,
MATPDE=4
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
void SUBR(create_matrix)(TYPE(crsMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  problem_t mat_type=FROM_FILE; // default: assume that 'which' is a file name

  int L; // problem size for Graphene (L x L grid) or the Anderson model (L^3 grid)

  // check if the first arg is of the grapheneN or andersonN form,
  // otherwise try to read a file
  mat_type=FROM_FILE;
  int pos=0;
  if( str_starts_with(problem,"spinSZ") )
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

  if (mat_type!=FROM_FILE)
  {
    // make sure all remaining characters are
    for (int i=pos;i<strlen(problem);i++)
    {
      if (problem[i]<'0' || problem[i]>'9')
      {
        mat_type=FROM_FILE;
        break;
      }
    }
  }

  if (mat_type!=FROM_FILE)
  {
    const char* strL=problem+pos;
    L=atoi(strL);
    if (L<0) mat_type=FROM_FILE;
  }
  
  if (mat_type==GRAPHENE)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: Graphene %d x %d\n",L,L);
    double gamma=0.2;
    ghost_lidx_t WL[2];
    WL[0] = L;
    WL[1] = L;
  
    ghost_gidx_t DIM = WL[0]*WL[1];
    matfuncs_info_t info;
    // set problem size
    crsGraphene( -2, WL, NULL, NULL);
    // set disorder
    crsGraphene( -5, NULL, NULL,&gamma);
    // get matrix info
    crsGraphene( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &crsGraphene, ierr), *ierr);
  }
  else if (mat_type==ANDERSON)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: Anderson %d x %d x %d\n",L,L,L);
    ghost_lidx_t LL=L;
  
    matfuncs_info_t info;
    anderson( -2, &LL, NULL, NULL);
    anderson( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &anderson, ierr), *ierr);
  }
  else if (mat_type==SPINSZ)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: spinSZ[%d]\n",L);


    ghost_gidx_t DIM;
    ghost_lidx_t conf_spinZ[3] = {L,L/2,0};
    SpinChainSZ( -2, conf_spinZ, &DIM, NULL);

    matfuncs_info_t info;
    SpinChainSZ( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc)(mat,comm,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &SpinChainSZ, ierr), *ierr);
  }
  else if (mat_type==MATPDE)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: MATPDE %d x %d\n", L, L);
    gidx_t nrows, ncols;
    lidx_t row_nnz;
    MATPDE_initDimensions(L, L, &nrows, &row_nnz);
    ncols = nrows;
    PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc)(mat, comm, 
          nrows, ncols, row_nnz, &MATPDE_rowFunc, ierr), *ierr);
  }
  else
  {
    PHIST_SOUT(PHIST_INFO,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(crsMat_read)(mat,comm,(char*)problem,ierr),*ierr);
  }
  
  return;
}
#else
void SUBR(create_matrix)(TYPE(crsMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* ierr)
{
  *ierr=-99;
}
#endif
