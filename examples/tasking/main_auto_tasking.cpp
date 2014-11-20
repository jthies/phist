#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>


#ifdef PHIST_KERNEL_LIB_GHOST
#include "phist_macros.h"
#include "phist_ghost_macros.hpp"
#include "phist_kernels.h"

void printHello(int n)
{
PHIST_GHOST_TASK_BEGIN
  std::cout << "Hello, World! ["<<n<<"]"<<std::endl;
PHIST_GHOST_TASK_END
}

int main(int argc, char* argv[])
{
  int ierr;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

  for (int i=1;i<=3;i++)
    {
    printHello(i);
    }

  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
}

#else
int main()
{
  std::cout << "this driver does nothing unless you use GHOST as kernel lib\n";
  return 0;
}
#endif
