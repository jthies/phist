#include "Checkpoint.hpp"
#include "UserClass.hpp"

int main()
{
  UserClass *my_object = new UserClass();
  int my_int=42;
  double my_double=42.9;
  Checkpoint backup;
  backup.add("i",&my_int);
  backup.add("x",&my_double);
  backup.add("my_object",my_object);
  backup.write();
}
