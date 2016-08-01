namespace phist
{

  namespace testing
  {

void SUBR(read_mat)(const char* filebase,phist_const_comm_ptr comm, int nglob,TYPE(sparseMat_ptr) *ptr, int* iflag)
{
  int iflag_in=*iflag;
  iflag_in |= PHIST_OUTLEV>=PHIST_DEBUG ? 0 : PHIST_SPARSEMAT_QUIET;
  *ptr = NULL;
  char tpc = ::phist::ScalarTraits< _ST_ >::type_char();
  char mmfile[256],hbfile[256],binfile[256],crsfile[256];
  sprintf(mmfile,"%c%s%d.mm",tpc,filebase,nglob);
  if (::phist::ScalarTraits< _ST_ >::is_complex())
  {
    sprintf(hbfile,"%c%s%d.cua",tpc,filebase,nglob);
  } 
  else
  {
    sprintf(hbfile,"%c%s%d.rua",tpc,filebase,nglob);
  }

  sprintf(binfile,"%c%s%d.bin",tpc,filebase,nglob);
  sprintf(crsfile,"%c%s%d.crs",tpc,filebase,nglob);
  
  PHIST_SOUT(PHIST_DEBUG, "Looking for matrix \'%s\'..\n", filebase);

  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", binfile);
  *iflag = iflag_in;
  SUBR(sparseMat_read_bin)(ptr,comm,binfile,iflag);
  if (*iflag==PHIST_SUCCESS) return;

  // try same format, different extension
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", crsfile);
  *iflag = PHIST_OUTLEV>=PHIST_DEBUG ? 0 : PHIST_SPARSEMAT_QUIET;
  SUBR(sparseMat_read_bin)(ptr,comm,crsfile,iflag);
  if (*iflag==PHIST_SUCCESS) return;
  
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", mmfile);
  *iflag = iflag_in;
  SUBR(sparseMat_read_mm)(ptr,comm,mmfile,iflag);
  if (*iflag==PHIST_SUCCESS) return;

  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", hbfile);
  *iflag = iflag_in;
  SUBR(sparseMat_read_hb)(ptr,comm,hbfile,iflag);
}//read_mat

// some functions for initializing matrices and mvecs
int PHIST_TG_PREFIX(idfunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
  *len=1;
  _ST_* val = (_ST_*)vval;
  val[0]=(_ST_)1.0;
  cols[0]=row;
  return 0;
}
/*
int PHIST_TG_PREFIX(some_rowFunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  _ST_* val = (_ST_*)vval;

  *len=5;
  for (int i=0; i<*len; i++)
  {
    cols[i]=(ghost_gidx)(((row+i-2)*3)%_N_);
    if (cols[i]<0) cols[i]+=_N_;
    val[i]=(ST)(i+1)/(ST)(row+1) + st::cmplx_I()*(ST)(row-cols[i]);
  }
  return 0;
}
*/
  int PHIST_TG_PREFIX(mvec123func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg)
  {
    _ST_* v= (_ST_*)val;
    int *int_arg=(int*)last_arg;
    int N  = int_arg[0];
    int NV = int_arg[1];
    *v = (_ST_)(i+1 + N*j);
    return 0;
  }

  int PHIST_TG_PREFIX(mvec321func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg)
  {
    _ST_* v= (_ST_*)val;
    int *int_arg=(int*)last_arg;
    int N  = int_arg[0];
    int NV = int_arg[1];
    *v = (_ST_)((N-i) + N*(NV-(j+1)));
    return 0;
  }


  }//testing
}//phist
