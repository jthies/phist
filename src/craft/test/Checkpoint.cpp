#include "Checkpoint.hpp"
#include "CpPOD.hpp"

void Checkpoint::write()
{
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
  {
    std::cout << "'"<<it->first << "' says hello like this: "<<std::endl;
    it->second->hello();
  }
}

