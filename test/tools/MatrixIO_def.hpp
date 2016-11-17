namespace phist
{

  namespace testing
  {

void SUBR(read_mat)(const char* filebase,phist_const_comm_ptr comm, 
        int nglob, int mglob, TYPE(sparseMat_ptr) *ptr, int* iflag)
{
  int iflag_in=*iflag;
  std::stringstream dimstr;
  dimstr<<nglob;
  if (mglob!=nglob) dimstr<<"x"<<mglob;
  iflag_in |= PHIST_OUTLEV>=PHIST_DEBUG ? 0 : PHIST_SPARSEMAT_QUIET;
  *ptr = NULL;
  char tpc = ::phist::ScalarTraits< _ST_ >::type_char();
  char mmfile[256],hbfile[256],binfile[256],crsfile[256];
  sprintf(mmfile,"%c%s%d.mm",tpc,filebase,nglob);
  if (::phist::ScalarTraits< _ST_ >::is_complex())
  {
    sprintf(hbfile,"%c%s%s.cua",tpc,filebase,dimstr.str().c_str());
  } 
  else
  {
    sprintf(hbfile,"%c%s%s.rua",tpc,filebase,dimstr.str().c_str());
  }

  sprintf(binfile,"%c%s%s.bin",tpc,filebase,dimstr.str().c_str());
  sprintf(crsfile,"%c%s%s.crs",tpc,filebase,dimstr.str().c_str());
  
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
  static phist_gidx gnrows=-1;
  static phist_gidx gncols=-1;
  if (row==-1)
  {
    gnrows=cols[0];
    gncols=cols[1];
    return 0;
  }
  if (gnrows>0 && row>=gnrows) return -1;
  if (row<0) return -1;
  if (gncols>0 && row>=gncols)
  {
    *len=0;
    return 0;
  }
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

int PHIST_TG_PREFIX(hpd_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg)
{
#include "phist_std_typedefs.hpp"
  static ghost_gidx gnrows=-1;
  _ST_ *vals=(_ST_*)vval;
  
  if (vals) vals[0]=st::one();
  if (cols && row>=0) cols[0]=row;

// create identity matrix (just while developing stuff with B-inner products
//  if (len) *len=1;
//  return 0;

  if (row<0)
  {
    gnrows=cols[0];
    return 0;
  }
  else if (gnrows<0)
  {
    PHIST_SOUT(PHIST_ERROR,"%s not correctly initialized, call with row=-1 and cols[0]=gnrows first!",__FUNCTION__);
    return -1;
  }
  else if (row==0)
  {
    *len=2;
    cols[1]=row+1;
#ifdef IS_COMPLEX
    vals[1]=(_ST_)-0.5*st::cmplx_I();
#else
    vals[1]=-0.5*st::one();
#endif
  }
  else if (row==gnrows-1)
  {
    *len=2;
    cols[1]=row-1;
#ifdef IS_COMPLEX
    vals[1]=(_ST_)+0.5*st::cmplx_I();
#else
    vals[1]=-0.5*st::one();
#endif
  }
  else
  {
    *len=3;
    cols[1]=row-1;
    cols[2]=row+1;
#ifdef IS_COMPLEX
    vals[1]=(_ST_)+0.5*st::cmplx_I();
    vals[2]=(_ST_)-0.5*st::cmplx_I();
#else
    vals[1]=-0.5*st::one();
    vals[2]=-0.5*st::one();
#endif
  }
  return 0;
}

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
