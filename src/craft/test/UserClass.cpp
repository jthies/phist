#include <iostream>
#include "UserClass.hpp"

UserClass::UserClass()
{
}

UserClass::~UserClass()
{
}

void UserClass::hello() const
{
  std::cout << "This is the user-defined class"<<std::endl;
}
