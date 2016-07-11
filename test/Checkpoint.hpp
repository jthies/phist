#ifndef CHECKPOINT_HPP
#define CHECKPOINT_HPP

#include <string>

#include "Checkpointable.hpp"
#include "CpPOD.hpp"

#include <string>
#include <map>

class Checkpoint
{
  public:
  
  typedef std::map<std::string,Checkpointable const*> cp_const_map;
  typedef std::map<std::string,Checkpointable *> cp_copy_map;
  
  void write();
  
  // contains points to user workspace
  cp_const_map objects;
  // contains asynchronous copies of objects
  cp_copy_map async_copies;

// implementation of add() for anything that is Checkpointable,
// anything else will give an error message
void add(std::string label, Checkpointable const* p)
{
    objects[label] = p;
}

// specialized implementation of add for POD (plain old data), which is obviously
// not derived from our base class Checkpointable. The same can be done for e.g. 
// ghost_densemat
// TODO: here a new object is created that will not be deleted unless we have some
//       memory management in the background, the "normal" add function just adds 
//       the pointer to the map and the person who created the object is responsible
//       for deleting it.
  void add(std::string label, int const* i){this->add(label, new CpPOD<int>(i));}
//  void add<int>(std::string label, int const* i){this->add(label, new CpPOD<int>(i));}

  void add(std::string label, double const* d){this->add(label, new CpPOD<double>(d));}

#endif

};
