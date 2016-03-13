#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <unistd.h>


#include "phist_macros.h"
#include "phist_tasks.h"
#include "phist_kernels.h"


void printHello(int n, int* iflag)
{
PHIST_TASK_DECLARE(printHelloTask)
PHIST_TASK_BEGIN(printHelloTask)
  std::cout << "Hello, World! ["<<n<<"]"<<std::endl;
PHIST_TASK_END(iflag)
}

void asyncPrintHello(int* iflag)
{
  int *iflag_t2 = new int;
  phist_map_ptr map = NULL;
  phist_comm_ptr comm = NULL;
  phist_Dmvec_ptr mvec = NULL;
  PHIST_CHK_IERR(phist_comm_create(&comm,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_create(&map,comm,(phist_gidx)1000,iflag),*iflag);
  PHIST_CHK_IERR(phist_Dmvec_create(&mvec,map,1,iflag),*iflag);

PHIST_TASK_DECLARE(Task1)
PHIST_TASK_DECLARE(Task2)
PHIST_TASK_DECLARE(Task3)

PHIST_TASK_BEGIN(Task1)
  sleep(1);
  std::cout << "Hello from task 1" << std::endl;
PHIST_TASK_END_NOWAIT(iflag);

PHIST_TASK_BEGIN(Task2)
  {
    //int *iflag_t2 = new int;
    *iflag_t2 = 0;
    PHIST_CHK_IERR(phist_Dmvec_put_value(mvec,1.,iflag_t2),*iflag_t2);
    double norm;
    PHIST_CHK_IERR(phist_Dmvec_norm2(mvec,&norm,iflag_t2),*iflag_t2);
    std::cout << "Hello from task 2, trying to produce an error!" << std::endl;
    PHIST_CHK_IERR(*iflag_t2=1,*iflag_t2);
  }
PHIST_TASK_END_NOWAIT(iflag)

PHIST_TASK_BEGIN(Task3)
  std::cout << "Hello from task 3" << std::endl;
PHIST_TASK_END_NOWAIT(iflag);

PHIST_TASK_WAIT(Task3,iflag)
PHIST_TASK_WAIT(Task2,iflag)
PHIST_TASK_WAIT(Task1,iflag)

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

