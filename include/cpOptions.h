#ifndef __CPOPTIONS_H__
#define __CPOPTIONS_H__

#include <string>

class CpOptions{
private:
  int cpFreq;
  int nIter;

public:
  CpOptions();
  ~CpOptions();

  void setCpFreq ( const int cpFreq_);
  int getCpFreq();
  void setnIter(const int nIter_);
  int getnIter();
  
};

#endif
