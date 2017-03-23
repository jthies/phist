#ifndef __ADDEDCPTYPES_HPP__
#define __ADDEDCPTYPES_HPP__

class Checkpoint;
#include "checkpoint.hpp"
#include "cpPOD.hpp"
#include "cpArray.hpp"
#ifdef GHOST_CP 
#include "cpTypes/cpGhost/cpGhost.hpp"
#endif

#ifdef PHIST_CP 
#include "cpTypes/cpPhistMvec/cpPhistMvec.cpp"
#include "cpTypes/cpPhistSdMat/cpPhistSdMat.cpp"
#endif

#ifdef MKL_CP
#include "mkl.h"
#include "cpTypes/cpMkl/cpMkl.hpp"
#include "cpTypes/cpMkl/cpMklArray.hpp"
#endif

// specialized implementation of add functions
// NOTE: here a new object is created that will not be deleted unless we have some
//       memory management in the background, the "normal" add function just adds 
//       the pointer to the map and the person who created the object is responsible
//       for deleting it.

// ===== POD ===== //
template <class T>
int addCpType(Checkpoint * cp, std::string label, T * const i){
  cp->addToMap(label, new CpPOD<T>(i, cp->getCpComm()));
}


// ===== POD ARRAY ===== // 
//TODO: note: if an array is defined statically i.e. int a[5]; the array-argument in add function should be preceded with (int *) a
template <class T>
int addCpType(Checkpoint * cp, std::string label, T* const arrayPtr_, const size_t nRows_){     
  cp->addToMap(label, new CpArray<T>(arrayPtr_, nRows_, cp->getCpComm()));
}

// ===== POD MULTI-ARRAY ===== // 
template <class T>
int addCpType(Checkpoint * cp, std::string label, T** const arrayPtr_, const size_t nRows_, const size_t nCols_, const int toCpCol_){
  cp->addToMap(label, new CpMultiArray<T>(arrayPtr_, nRows_, nCols_, toCpCol_, cp->getCpComm() ));
}

#ifdef GHOST_CP
// ===== GHOST DENSE MATRIX ===== //
int addCpType(Checkpoint * cp, std::string label, ghost_densemat * const GDM)
{
  cp->addToMap(label, new CpGhostDenseMat(GDM, cp->getCpComm() ));
}

int addCpType(Checkpoint * cp, std::string label, ghost_densemat ** const GDMArray, const size_t nDenseMat_, const int toCpDenseMat_)
{
  cp->addToMap(label, new CpGhostDenseMatArray(GDMArray, nDenseMat_, toCpDenseMat_, cp->getCpComm() ) );
}

// ===== GHOST SPARSE MATRIX ===== // TODO: add this functionality if needed by users
//int Checkpoint::add(std::string label, ghost_sparsemat * const GSM)
//{
//		this->add(label, new CpGhostSparseMat(GSM));
//}

#endif

#ifdef PHIST_CP 
// ===== PHIST MVEC ===== // 
int addCpType(Checkpoint * cp, std::string label, TYPE(mvec_ptr) const mvec)
{	
  cp->addToMap(label, new TYPE(CpPhistMvec)(mvec) );
}

// ===== PHIST SDMAT ===== //	TODO: as MVEC, and SDMAT are both void*, they should be differenciated in some better way.  
int addCpType(Checkpoint * cp, std::string label, TYPE(sdMat_ptr) const sdMat, TYPE(sdMat_ptr) const temp)
{	
  cp->addToMap(label, new TYPE(CpPhistSdMat)(sdMat, cp->getCpComm()) );
}
#endif


#ifdef MKL_CP
// ===== MKL_Complex8 ===== // 
int addCpType(Checkpoint * cp, std::string label, MKL_Complex8 * const dataPtr){
  cp->addToMap(label, new CpMklComplex8(dataPtr, cp->getCpComm()));
}

// ===== MKL_Complex16 ===== // 
int addCpType(Checkpoint * cp, std::string label, MKL_Complex16 * const dataPtr){
  cp->addToMap(label, new CpMklComplex16(dataPtr, cp->getCpComm()));
}

// ===== MKL_Complex8 * (Array) ===== // 
int addCpType(Checkpoint * cp, std::string label, MKL_Complex8 * const dataPtr, const size_t nElem){
  cp->addToMap(label, new CpMklComplex8Array(dataPtr, nElem, cp->getCpComm()));
}

// ===== MKL_Complex16 * (Array) ===== // 
int addCpType(Checkpoint * cp, std::string label, MKL_Complex16 * const dataPtr, const size_t nElem){
  cp->addToMap(label, new CpMklComplex16Array(dataPtr, nElem, cp->getCpComm()));
}
#endif


#endif

