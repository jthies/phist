#include <cstdlib>
#include <fstream>
#include <iostream>


#include "include/cpTypes/cpGhost/cpGhost.hpp"
#include "cpHelperFuncs.hpp"

CpGhostDenseMat::CpGhostDenseMat(ghost_densemat *  dataPtr_, const MPI_Comm cpMpiComm_){
		dataPtr = dataPtr_;
		ghost_densemat_create(&asynData, dataPtr->map, dataPtr->traits); 
		ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);


}

int CpGhostDenseMat::update(){
  ghost_densemat_init_densemat(asynData, dataPtr, 0, 0);
  return EXIT_SUCCESS;
}

int CpGhostDenseMat::write( const std::string * filename){
  ghost_densemat_to_file(asynData,(char *) (*filename).c_str(), MPI_COMM_SELF);
  return EXIT_SUCCESS;
}

int CpGhostDenseMat::read(const std::string * filename){
  ghost_densemat_init_file(dataPtr,(char *) (*filename).c_str(), MPI_COMM_SELF);
  return EXIT_SUCCESS;
}



CpGhostDenseMatArray::CpGhostDenseMatArray(ghost_densemat **  dataPtr_, const size_t nDenseMat_, const int toCpDenseMat_, const MPI_Comm cpMpiComm_){
  dataPtr = dataPtr_;
  nDenseMat = nDenseMat_;
  toCpDenseMat = toCpDenseMat_;
  cyclicCpCounter = 0;
  
  asynData = new ghost_densemat*[nDenseMat];
  for(size_t i = 0; i < nDenseMat ; ++i)
  {
    ghost_densemat_create( &asynData[i], dataPtr[i]->map, dataPtr[i]->traits);
    ghost_densemat_init_densemat ( asynData[i], dataPtr[i], 0, 0);	
  }
}


int CpGhostDenseMatArray::update(){
  for(size_t i = 0; i < nDenseMat ; ++i)
  {
    ghost_densemat_init_densemat ( asynData[i], dataPtr[i], 0, 0);	
  }
  return EXIT_SUCCESS;
}

int CpGhostDenseMatArray::write( const std::string * filename){
  if(toCpDenseMat == ALL){
    craftDbg(3, "CpGhostDenseMatArray::write() toCpDenseMat is ALL");
    for(size_t i = 0; i < nDenseMat ; ++i)
    {
      ghost_densemat_to_file(asynData[i],(char *) (*filename).c_str(), MPI_COMM_SELF);
    }
  }
  else if(toCpDenseMat == CYCLIC){		// TODO: testing to be done
    craftDbg(3, "CpGhostDenseMatArray::write() cyclicCpCounter = %d", cyclicCpCounter);
    ghost_densemat_to_file(asynData[cyclicCpCounter],(char *) (*filename).c_str(), MPI_COMM_SELF);
    // ===== write the metadata file for cyclicCpCounter ===== // 
    std::string filenameMD;
    filenameMD = *filename + ".metadata";
    std::ofstream fstrMD;
    fstrMD.open ((filenameMD).c_str(), std::ios::out | std::ios::binary );	
    if(fstrMD.is_open()){
      fstrMD.write( (char *)&(cyclicCpCounter), sizeof (size_t) );
      fstrMD.close();
    }
    ++cyclicCpCounter;
    if( cyclicCpCounter == nDenseMat ){
      cyclicCpCounter = 0;	
    }
  }
  else if(toCpDenseMat	>= 0){
    craftDbg(3, "CpGhostDenseMatArray::write() toCpDenseMat is: %d", toCpDenseMat);
    ghost_densemat_to_file(asynData[toCpDenseMat],(char *) (*filename).c_str(), MPI_COMM_SELF);
  }
  else{
    std::cerr << "ERROR at " << __FILE__ << " " << __LINE__ << std::endl;			
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int CpGhostDenseMatArray::read(const std::string * filename){
  if(toCpDenseMat == ALL){
    craftDbg(3, "CpGhostDenseMatArray::read(): toCpDenseMat is ALL");
    for(size_t i = 0; i < nDenseMat ; ++i)
    {
      ghost_densemat_init_file(dataPtr[i],(char *) (*filename).c_str(), MPI_COMM_SELF);
    }
  }
  else if(toCpDenseMat == CYCLIC){			// TODO: testing to be done
      // ===== read the metadata file for cyclicCpCounter ===== // 
      std::string filenameMD;
      filenameMD = *filename + ".metadata";
      std::ifstream fstrMD;
      fstrMD.open ((filenameMD).c_str(), std::ios::in | std::ios::binary );	
      if(fstrMD.is_open()){
        fstrMD.read( (char *)&(cyclicCpCounter), sizeof (size_t) );
        fstrMD.close();
      }
      craftDbg(3, "CpGhostDenseMatArray::read(): cyclicCpCounter is = %d", cyclicCpCounter);
      ghost_densemat_init_file(dataPtr[cyclicCpCounter],(char *) (*filename).c_str(), MPI_COMM_SELF);
      ++cyclicCpCounter;
      if( cyclicCpCounter == nDenseMat ){
        cyclicCpCounter = 0;	
      }
  }
  else if( toCpDenseMat >= 0 ){
    craftDbg(3, "CpGhostDenseMatArray::read(): toCpDenseMat is = %d", toCpDenseMat);
    ghost_densemat_init_file(dataPtr[toCpDenseMat],(char *) (*filename).c_str(), MPI_COMM_SELF);
  }
  else{
    std::cerr << "ERROR at " << __FILE__ << " " << __LINE__ << std::endl;			
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

CpGhostSparseMat::CpGhostSparseMat(ghost_sparsemat *  dataPtr_){
  dataPtr = dataPtr_;
  std::cerr << "CpGhostSparseMat checkpoints are not supported yet." << std::endl;
}

int CpGhostSparseMat::update(){
  std::cerr << "CpGhostSparseMat checkpoints are not supported yet." << std::endl;
  return EXIT_FAILURE;
}

int CpGhostSparseMat::write( const std::string * filename){
  std::cerr << "CpGhostSparseMat checkpoints are not supported yet." << std::endl;
  return EXIT_FAILURE;
}

int CpGhostSparseMat::read(const std::string * filename){
  std::cerr << "CpGhostSparseMat checkpoints are not supported yet." << std::endl;
  return EXIT_FAILURE;
}



