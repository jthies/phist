#include <cstdint>
#include "phist_config.h"
#include "phist_macros.h"

#ifdef PHIST_HAVE_COLPACK
#include "ColPack/ColPackHeaders.h"
#endif
// interface to ColPack's graph coloring algorithms. 
// Input: local CRS arrays as used in the builtin kernels (only local edges are considered)
//        dist: only distance 2 supported for now
//        *ncolors: will contain the number of colurs
//        colors: dimension nrows, preallocated by user, contains 0-based color values.
extern "C" void colpack_v1_0_8(int nrows, int64_t* row_ptr, int64_t* nonlocal_ptr, int* col_idx, 
        int* ncolors, int* colors, int dist, int idx_base, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#ifndef PHIST_HAVE_COLPACK
  PHIST_SOUT(PHIST_ERROR,"%s not implemented, ColPack not available.\n"
                         "(file %s, line %d)\n",
  __FUNCTION__,__FILE__,__LINE__);
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
#else
  *iflag=0;

  // ColPack uses a single long sequence of inherited classes, e.g.
  // GraphColoring 
  //      GraphOrdering 
  //              GraphInputOutput 
  //              ...
  // ... all the way, so we just need a single object to represent the graph and the ordering 
  // and everything
  ColPack::GraphColoring *GC=new ColPack::GraphColoring();

  // the input format this package accepts is called "ADOL-C compressed row". I have
  // looked in an obscure paper, ADOL-C is something about automatic differentiation
  // and tapes, and their format is this:
  uint32_t** adolc = new uint32_t*[nrows]; 
                          // each of these pointers points to row i, but 
                          // entry 0 is the length of row i, so we need to
                          // re-allocate and copy instead of passing in CRS

 // first we need to count the local nonzeros (TODO - or are they readily 
 // available in crsmat_module?)
 int64_t nzloc=0;
 for (int i=0;i<nrows;i++) 
 {
   nzloc+=nonlocal_ptr[i]-row_ptr[i];
 }

 uint32_t *adolc_data=new uint32_t[nzloc+nrows];

//#pragma omp parallel for schedule(static)
 int64_t pos=0;
 for (int i=0;i<nrows;i++)
 {
   adolc[i]=&(adolc_data[pos]);
   // first entry in row is the local row length
   adolc_data[pos++]=nonlocal_ptr[i]-row_ptr[i];
   for (int j=row_ptr[i];j<nonlocal_ptr[i];j++)
   {
     // subtract idx_base to account for 1-based input
     adolc_data[pos++]=col_idx[j-idx_base]-idx_base;
   }
 }

  int i_HighestDegree = GC->BuildGraphFromRowCompressedFormat(adolc, nrows);
  
  delete [] adolc_data;
  delete [] adolc;

  if (dist==1)
  {
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  }
  if (dist==2)
  {
    // this function returns 1 on success
    PHIST_CHK_IERR(*iflag=GC->DistanceTwoColoring()-1,*iflag);
#if PHIST_OUTLEV>=PHIST_VERBOSE
int verbose=2;
#else
int verbose=1;
#endif
#ifdef TESTING__disabled_because_of_issue_67
    PHIST_SOUT(PHIST_VERBOSE,"thoroughly test local dist-2 coloring\n");
    PHIST_CHK_IERR(*iflag=GC->CheckDistanceTwoColoring(verbose),*iflag);
#else
    PHIST_SOUT(PHIST_VERBOSE,"quickly test local dist-2 coloring\n");
    PHIST_CHK_IERR(*iflag=GC->CheckQuickDistanceTwoColoring(verbose),*iflag);
#endif
  }
  else
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }

  *ncolors = GC->GetVertexColorCount();
  std::vector<int>& colorVec=*(GC->GetVertexColorsPtr());

  PHIST_CHK_IERR(*iflag=nrows-colorVec.size(),*iflag);
  
  for (int i=0;i<nrows;i++)
  {
    colors[i]=colorVec[i];
  }

PHIST_SOUT(PHIST_VERBOSE,"number of local colors (rank 0): %d\n",*ncolors);

  delete GC;
  return;
#endif
}
