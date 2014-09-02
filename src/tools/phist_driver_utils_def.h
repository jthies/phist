// auto-detects the file type by looking at the file extension
void SUBR(crsMat_read)(TYPE(crsMat_ptr)* A, char* filename, int* ierr)
{
  char *c=&filename[0];
  int isMM=0;
  int isCRS=0;
  int isHB=0;
  
  while (*c)
  {
    if (*c=='.') break;
    c++;
  }
  c++;
  isMM=!strcmp(c,"mm");
  isCRS=!strcmp(c,"crs");
  isCRS=isCRS||!strcmp(c,"bin");
  isHB=!strcmp(c,"rua")||!strcmp(c,"cua");
  if (isMM)
  {
    PHIST_CHK_IERR(SUBR(crsMat_read_mm)(A,filename,ierr),*ierr);
  }
  else if (isHB)
  {
    PHIST_CHK_IERR(SUBR(crsMat_read_hb)(A,filename,ierr),*ierr);
  }
  else if (isCRS)
  {
    PHIST_CHK_IERR(SUBR(crsMat_read_bin)(A,filename,ierr),*ierr);
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

#include "matfuncs.h"

typedef enum {
FROM_FILE=0,
GRAPHENE=1,
ANDERSON=2
} problem_t;

void SUBR(create_matrix)(TYPE(crsMat_ptr)* mat, const char* problem, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  problem_t mat_type=FROM_FILE; // default: assume that 'which' is a file name

  int L; // problem size for Graphene (L x L grid) or the Anderson model (L^3 grid)

  // check if the first arg is of the grapheneN or andersonN form,
  // otherwise try to read a file
  mat_type=FROM_FILE;

  if (strlen(problem)>=8)
  {
    if (strncmp(problem,"graphene",8)==0)
    {
      mat_type=GRAPHENE;
    }
    else if (strncmp(problem,"anderson",8)==0)
    {
      mat_type=ANDERSON;
    }
    if (mat_type!=FROM_FILE)
    {
      // make sure all remaining characters are
      for (int i=8;i<strlen(problem);i++)
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
      const char* strL=problem+8;
      L=atoi(strL);
      if (L<0) mat_type=FROM_FILE;
    }
  }
  
  if (mat_type==GRAPHENE)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: Graphene %d x %d\n",L,L);
    ghost_lidx_t WL[2];
    WL[0] = L;
    WL[1] = L;
  
    ghost_gidx_t DIM = WL[0]*WL[1];
    matfuncs_info_t info;
    crsGraphene( -2, WL, NULL, NULL);
    crsGraphene( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc)(mat,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &crsGraphene, ierr), *ierr);
  }
  else if (mat_type==ANDERSON)
  {
    PHIST_SOUT(PHIST_INFO,"problem type: Anderson %d x %d %d\n",L,L,L);
    ghost_lidx_t LL=L;
  
    matfuncs_info_t info;
    anderson( -2, &LL, NULL, NULL);
    anderson( -1, NULL, NULL, &info);

    PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc)(mat,
        (gidx_t)info.nrows, (gidx_t)info.ncols, (lidx_t)info.row_nnz,
        &anderson, ierr), *ierr);
  }
  else
  {
    PHIST_SOUT(PHIST_INFO,"read matrix from file '%s'\n",problem);
    PHIST_CHK_IERR(SUBR(crsMat_read)(mat,(char*)problem,ierr),*ierr);
  }
  
  return;
}
#else
void SUBR(create_matrix)(TYPE(crsMat_ptr)* mat, const char* problem, int* ierr)
{
  *ierr=-99;
}
#endif
