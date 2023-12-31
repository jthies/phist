/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include <cstdint>
#include "phist_config.h"
#include "phist_macros.h"
#ifndef PHIST_TESTING
#define PHIST_TESTING
#endif
#ifdef PHIST_HAVE_COLPACK
#include "ColPackHeaders.h"
#endif
// interface to ColPack's graph coloring algorithms. 
// Input: CRS arrays as used in the builtin kernels (off-processor elements are renumbered and bipartite row-coloring 
//        is used so that all local rows of a color have a distance at least of 'dist' taking ghost nodes into account)
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
  ColPack::BipartiteGraphPartialColoring *GC=new ColPack::BipartiteGraphPartialColoring();

  // the input format this package accepts is called "ADOL-C compressed row". I have
  // looked in an obscure paper, ADOL-C is something about automatic differentiation
  // and tapes, and their format is this:
  uint32_t** adolc = new uint32_t*[nrows]; 
                          // each of these pointers points to row i, but 
                          // entry 0 is the length of row i, so we need to
                          // re-allocate and copy instead of passing in CRS

 // first we need to count the local nonzeros and determine the number of columns
 // in the local part of the matrix. If there are nrows local rows and the largest
 // off-processor column entry is c_max, ncols=nrows+c_max. Since this coloring is
 // only local to the MPI process, we do not need a consistent numbering across
 // MPI processes.
 int64_t nzloc=0;
 int ncols=nrows;
 int nthreads=1;
#pragma omp parallel
{
#pragma omp master
 nthreads=omp_get_num_threads();
}

 int64_t thread_nrows[nthreads];
 int64_t thread_nzloc[nthreads];
 int64_t thread_c_max[nthreads];
 
 for (int i=0; i<nthreads;i++)
 {
   thread_nrows[i]=0;
   thread_nzloc[i]=0;
   thread_c_max[i]=0;
 }
 
#pragma omp parallel
{
  int thread_id=omp_get_thread_num();
#pragma omp for schedule(static)
  for (int i=0;i<nrows;i++)
  {
    thread_nrows[thread_id]++;
    thread_nzloc[thread_id]+=row_ptr[i+1]-row_ptr[i];
    for (int j=nonlocal_ptr[i]; j<row_ptr[i+1]; j++)
    {
      thread_c_max[thread_id] = std::max(thread_c_max[thread_id], (int64_t)col_idx[j-idx_base]);
    }
  }
}

for (int i=0;i<nthreads;i++)
{
  nzloc+=thread_nzloc[i];
  ncols+=thread_c_max[i];
}


 uint32_t *adolc_data=new uint32_t[nzloc+nrows];

#pragma omp parallel
{
  int thread_id=omp_get_thread_num();
  int64_t pos=0;
  for (int i=0; i<thread_id; i++)
  {
    pos+=thread_nrows[i]+thread_nzloc[i];
  }
#pragma omp for schedule(static)
  for (int i=0;i<nrows;i++)
  {
    adolc[i]=&(adolc_data[pos]);
    // first entry in row is the local row length
    adolc_data[pos++]=row_ptr[i+1]-row_ptr[i];
    for (int j=row_ptr[i];j<nonlocal_ptr[i];j++)
    {
      // subtract idx_base to account for 1-based input
      adolc_data[pos++]=col_idx[j-idx_base]-idx_base;
    }
    for (int j=nonlocal_ptr[i];j<row_ptr[i+1];j++)
    {
      // add nrows to assign unique node IDs to the off-processor column indices.
      adolc_data[pos++]=col_idx[j-idx_base]+nrows;
    }
  }
}

  int i_HighestDegree = GC->BuildBPGraphFromRowCompressedFormat(adolc, nrows,ncols);
  
  delete [] adolc_data;
  delete [] adolc;

  if (dist==1)
  {
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  }
  if (dist==2)
  {
    // this function returns 1 on success
    PHIST_CHK_IERR(*iflag=GC->PartialDistanceTwoRowColoring()-1,*iflag);
#if PHIST_OUTLEV>=PHIST_VERBOSE
int verbose=2;
#else
int verbose=1;
#endif
//#ifdef PHIST_TESTING
    PHIST_SOUT(PHIST_VERBOSE,"thoroughly test local dist-2 coloring\n");
    PHIST_CHK_IERR(*iflag=GC->CheckPartialDistanceTwoColoring()-1,*iflag);
/*
#else
    PHIST_SOUT(PHIST_VERBOSE,"quickly test local dist-2 coloring\n");
    PHIST_CHK_IERR(*iflag=GC->CheckQuickDistanceTwoColoring(verbose),*iflag);
#endif
*/
  }
  else
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }

  *ncolors = GC->GetVertexColorCount();
  std::vector<int> colorVec;
  GC->GetVertexPartialColors(colorVec);

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
