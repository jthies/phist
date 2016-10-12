#ifndef __CPOPTIONS_H__
#define __CPOPTIONS_H__

#include <string>

class CpOptions{
private:
  int cpFreq;
  int nIter;
  std::string cpPath;
  bool isRestarted;
  bool withSCR;

public:
  CpOptions();
  ~CpOptions();
  void setCpPath ( const std::string cpPath_); 

  void setCpFreq ( const int cpFreq_);
  int getCpFreq();
  void setnIter(const int nIter_);
  int getnIter();
  std::string getCpPath(); 
  
  void setRestartStatus(bool status);
  bool getRestartStatus();

};

#endif
