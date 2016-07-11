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
template <typename T>
void add(std::string label, T const* p)
{
  Checkpointable const* cp_p = dynamic_cast<Checkpointable const*>(p);
  if (cp_p==NULL) 
  {
                  std::cerr << "can't add object with label '"<<label<<"'\n"
                            << "because it is not derived from Checkpointable\n"
                            << "and has no specialized implementation.\n";
    return;
  }
  else
  {
    objects[label] = cp_p;
  }
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
