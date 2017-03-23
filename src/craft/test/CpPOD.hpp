#ifndef CP_POD_HPP
#define CP_POD_HPP

#include <iostream>

#include "Checkpointable.hpp"

template<typename T>
class CpPOD : public Checkpointable
{
  public:
  
  CpPOD(T const* x)  : x_(x) {}

  ~CpPOD(){}
  
  void hello() const
  {
    std::cout << "Hi, this is a POD with value "<<*x_<<std::endl;
  }

  protected:
  
  T const* x_;

};


#endif
