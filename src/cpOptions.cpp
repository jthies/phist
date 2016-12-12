#include <stdio.h>
#include "include/cpOptions.h"



CpOptions::CpOptions(){
  cpFreq = -9999999;
  nIter = 0;
}

CpOptions::~CpOptions(){
  
}

void CpOptions::setCpFreq(const int cpFreq_){
  cpFreq = cpFreq_;
}

void CpOptions::setnIter(const int nIter_){
  nIter = nIter_;
}

int CpOptions::getCpFreq (){
  return cpFreq;
}

int CpOptions::getnIter(){
  return nIter;
}

