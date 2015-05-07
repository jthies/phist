#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <unistd.h>


#ifdef PHIST_KERNEL_LIB_GHOST
#include "phist_macros.h"
#include "phist_ghost_macros.hpp"
#include "phist_kernels.h"


void printHello(int n, int* iflag)
{
PHIST_GHOST_TASK_BEGIN(printHelloTask)
  std::cout << "Hello, World! ["<<n<<"]"<<std::endl;
PHIST_GHOST_TASK_END(iflag)
}

void asyncPrintHello(int* iflag)
{
  int *iflag_t2 = new int;
  map_ptr_t map = NULL;
  comm_ptr_t comm = NULL;
  Dmvec_ptr_t mvec = NULL;
  PHIST_CHK_IERR(phist_comm_create(&comm,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_create(&map,comm,(gidx_t)1000,iflag),*iflag);
  PHIST_CHK_IERR(phist_Dmvec_create(&mvec,map,1,iflag),*iflag);

PHIST_GHOST_TASK_BEGIN(Task1)
  sleep(1);
  std::cout << "Hello from task 1" << std::endl;
PHIST_GHOST_TASK_END_NOWAIT(iflag);

PHIST_GHOST_TASK_BEGIN(Task2)
  {
    //int *iflag_t2 = new int;
    *iflag_t2 = 0;
    PHIST_CHK_IERR(phist_Dmvec_put_value(mvec,1.,iflag_t2),*iflag_t2);
    double norm;
    PHIST_CHK_IERR(phist_Dmvec_norm2(mvec,&norm,iflag_t2),*iflag_t2);
    std::cout << "Hello from task 2, trying to produce an error!" << std::endl;
    PHIST_CHK_IERR(*iflag_t2=1,*iflag_t2);
  }
PHIST_GHOST_TASK_END_NOWAIT(iflag)

PHIST_GHOST_TASK_BEGIN(Task3)
  std::cout << "Hello from task 3" << std::endl;
PHIST_GHOST_TASK_END_NOWAIT(iflag);

PHIST_GHOST_TASK_WAIT(Task3,iflag)
PHIST_GHOST_TASK_WAIT(Task2,iflag)
PHIST_GHOST_TASK_WAIT(Task1,iflag)

  // check the error of Task2
  //PHIST_CHK_IERR(*iflag = *iflag_t2-1, *iflag);

  PHIST_CHK_IERR(phist_Dmvec_delete(mvec,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_delete(map,iflag),*iflag);
  PHIST_CHK_IERR(phist_comm_delete(comm,iflag),*iflag);
  delete iflag_t2;
}

int main(int argc, char* argv[])
{
  int iflag;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  PHIST_SOUT(PHIST_INFO, "synchronous tasks\n");
  for (int i=1;i<=3;i++)
  {
    PHIST_ICHK_IERR(printHello(i,&iflag),iflag);
  }

  PHIST_SOUT(PHIST_INFO, "asynchronous tasks\n");
  PHIST_ICHK_IERR(asyncPrintHello(&iflag),iflag);

  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
}

#else
int main()
{
  std::cout << "this driver does nothing unless you use GHOST as kernel lib\n";
  return 0;
}
#endif
