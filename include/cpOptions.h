#ifndef __CPOPTIONS_H__
#define __CPOPTIONS_H__

#include <string>

class CpOptions{
private:
  int cpFreq;
  int nIter;
  bool isRestarted;
  bool withSCR;

public:
  CpOptions();
  ~CpOptions();

  void setCpFreq ( const int cpFreq_);
  int getCpFreq();
  void setnIter(const int nIter_);
  int getnIter();
  
  void setRestartStatus(bool status);
  bool getRestartStatus();

};

#endif
