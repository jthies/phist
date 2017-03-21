#ifndef USER_CLASS_HPP
#define USER_CLASS_HPP

#include "Checkpointable.hpp"

class UserClass : public Checkpointable
{
  public:
  
  UserClass();
  
  ~UserClass();
  
  void hello() const;
};

#endif

