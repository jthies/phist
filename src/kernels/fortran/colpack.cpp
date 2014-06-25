#include <cstdint>

#include "phist_config.h"
#include "phist_macros.h"

#ifdef PHIST_HAVE_COLPACK
#include "ColPack/ColPackHeaders.h"
#endif
// interface to ColPack's graph coloring algorithms. 
// Input: local CRS arrays as used in the fortran kernels (only local edges are considered)
//        dist: only distance 2 supported for now
//        *ncolors: will contain the number of colurs
//        colors: dimension nrows, preallocated by user, contains 0-based color values.
extern "C" void colpack_v1_0_8(int nrows, int64_t* row_ptr, int64_t* nonlocal_ptr, int* col_idx, 
        int* ncolors, int* colors, int dist, int idx_base, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#ifndef PHIST_HAVE_COLPACK
  PHIST_SOUT(PHIST_ERROR,"%s not implemented, ColPack not available.\n"
                         "(file %s, line %d)\n",
  __FUNCTION__,__FILE__,__LINE__);
  *ierr=PHIST_NOT_IMPLEMENTED;
  return;
#else
  *ierr=0;
  // ColPack uses a single long sequence of inherited classes, e.g.
  // GraphColoring 
  //      GraphOrdering 
  //              GraphInputOutput 
  //              ...
  // ... all the way, so we just need a single object to represent the graph and the ordering 
  // and everything
  ColPack::GraphColoring GC;

  // the input format this package accepts is called "ADOL-C compressed row". I have
  // looked in an obscure paper, ADOL-C is something about automatic differentiation
  // and tapes, and their format is this:
  uint32_t** adolc = new uint32_t*[nrows]; 
                          // each of these pointers points to row i, but 
                          // entry 0 is the length of row i, so we need to
                          // re-allocate and copy instead of passing in CRS

 // first we need to count the local nonzeros (TODO - or are they readily 
 // available in crsmat_module?)
 uint32_t nzloc=0;
 for (int i=0;i<nrows;i++) 
 {
   nzloc+=nonlocal_ptr[i]-row_ptr[i];
 }

 uint32_t *adolc_data=new uint32_t[nzloc+nrows];

//#pragma omp parallel for schedule(static)
 for (int i=0;i<nrows;i++)
 {
   // +i because we store the row length in pos 0 of each row
   int64_t pos=row_ptr[i]-row_ptr[0]+i; 
   adolc[i]=&(adolc_data[pos]);
   adolc_data[pos++]=nonlocal_ptr[i]-row_ptr[i];
   for (int j=row_ptr[i];j<nonlocal_ptr[i];j++)
   {
     adolc_data[pos++]=col_idx[j-row_ptr[0]]-idx_base;
   }
 }

  int i_HighestDegree = GC.BuildGraphFromRowCompressedFormat(adolc, nrows);
  
  delete [] adolc_data;
  delete [] adolc;

  if (dist==1)
  {
    PHIST_CHK_IERR(*ierr=PHIST_NOT_IMPLEMENTED,*ierr);
  }
  if (dist==2)
  {
    // this function returns 1 on success
    PHIST_CHK_IERR(*ierr=GC.DistanceTwoColoring()-1,*ierr);
#if PHIST_OUTLEV>=PHIST_VERBOSE
int verbose=2;
#else
int verbose=1;
#endif
#ifdef TESTING
    PHIST_CHK_IERR(*ierr=GC.CheckDistanceTwoColoring(verbose),*ierr);
#else
    PHIST_CHK_IERR(*ierr=GC.CheckQuickDistanceTwoColoring(verbose),*ierr);
#endif
  }
  else
  {
    PHIST_CHK_IERR(*ierr=PHIST_INVALID_INPUT,*ierr);
  }

  *ncolors = GC.GetVertexColorCount();
  std::vector<int>& colorVec=*(GC.GetVertexColorsPtr());

  PHIST_CHK_IERR(*ierr=nrows-colorVec.size(),*ierr);
  
  for (int i=0;i<nrows;i++)
  {
    colors[i]=colorVec[i];
  }
  return;
#endif
}
