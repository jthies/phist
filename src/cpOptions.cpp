#include "include/cpOptions.h"
#include <stdio.h>


CpOptions::CpOptions(){
  cpFreq = 0;
  nIter = 0;
  cpPath = new char[256];
  isRestarted = false;  
}


CpOptions::~CpOptions(){
  
}

void CpOptions::setCpFreq(const int cpFreq_){
  cpFreq = cpFreq_;
}

void CpOptions::setnIter(const int nIter_){
  nIter = nIter_;
}

void CpOptions::setCpPath(const std::string cpPath_){
  cpPath = cpPath_;
}

void CpOptions::setRestartStatus( bool status_){
  isRestarted = status_;
}

int CpOptions::getCpFreq (){
  return cpFreq;
}

std::string CpOptions::getCpPath (){
  return cpPath;
}

bool  CpOptions::getRestartStatus(){
  return isRestarted;
}

int CpOptions::getnIter(){
  return nIter;
}

