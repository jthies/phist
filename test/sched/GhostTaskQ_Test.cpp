#include <stdio.h>
#include <unistd.h>
#include <omp.h>

#include "ghost_taskq.h"
#include "ghost_util.h"

static void *task1(void *arg) {
#pragma omp parallel
#pragma omp single
	printf("task1: has %d threads\n",omp_get_num_threads());

	omp_set_num_threads(2);
#pragma omp parallel
	ghost_setCore(omp_get_thread_num());

#pragma omp parallel
	printf("task1: thread %d running @ core %d\n",omp_get_thread_num(), ghost_getCore());

	return NULL;
}


static void *task2(void *arg) {
#pragma omp parallel
#pragma omp single
	printf("task2: has %d threads\n",omp_get_num_threads());

	omp_set_num_threads(2);
#pragma omp parallel
	ghost_setCore(omp_get_thread_num()+2);

#pragma omp parallel
	printf("task2: thread %d running @ core %d\n",omp_get_thread_num(), ghost_getCore());

	return NULL;
}

static void *accuFunc(void *arg) 
{
        int *ret = (int *)ghost_malloc(sizeof(int));
        *ret = (*(int *)arg)+1;
        
#pragma omp parallel
        {
#pragma omp single
        printf("    ######### accuFunc: numThreads: %d\n",omp_get_num_threads());
        printf("    ######### accuFunc: thread %d running @ core %d\n",omp_get_thread_num(), ghost_getCore());
        }

        return ret;
}

static void *shortRunningFunc(void *arg) 
{

        usleep(5e5); // sleep 500 ms
        if (arg != NULL)
                printf("    ######### shortRunningFunc arg: %s\n",(char *)arg);

#pragma omp parallel
        {
#pragma omp single
        printf("    ######### shortRunningFunc: numThreads: %d\n",omp_get_num_threads());
        printf("    ######### shortRunningFunc: thread %d running @ core %d\n",omp_get_thread_num(), ghost_getCore());
        }

        return NULL;
}

static void *longRunningFunc(void *arg) 
{
        usleep(1e6); // sleep 1 sec
        if (arg != NULL)
                printf("    ######### longRunningFunc arg: %s\n",(char *)arg);

#pragma omp parallel
        {
#pragma omp single
        printf("    ######### longRunningFunc: numThreads: %d\n",omp_get_num_threads());
        printf("    ######### longRunningFunc: thread %d running @ core %d\n",omp_get_thread_num(), ghost_getCore());
        }

        return NULL;
}


class GhostTaskQ_Test: public SchedTest
  {
  public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    SchedTest::SetUp();
    
    t1 = (pthread_t *)ghost_malloc(sizeof(pthread_t));
    t2 = (pthread_t *)ghost_malloc(sizeof(pthread_t));
	
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    free(p1);
    free(p2);
    SchedTest::TearDown();
    }

  pthread_t* p1, p2;
  };

// this is a simple test from ghost itself, it doesn't really enqueue anything
TEST_F(GhostTaskQ_Test,OmpTest)
 {
    pthread_create(t1,NULL,task1,NULL);
    pthread_create(t2,NULL,task2,NULL);
 
    pthread_join(*t1,NULL);
    pthread_join(*t2,NULL);
}

TEST_F(GhostTaskQ_Test_WithQueue)
  {
        int foo = 42;

        printf("The thread pool consists of %d threads in %d locality domains\n",ghost_thpool->nThreads,ghost_thpool->nLDs);

        ghost_task_t *accuTask;
        ghost_task_t *lrTask;
        ghost_task_t *srTask;

        accuTask = ghost_task_init(1, 0, &accuFunc, &foo, GHOST_TASK_DEFAULT);

        printf("checking for correct task state and whether accuTask returns %d\n",foo+1);
        printf("state should be invalid: %s\n",ghost_task_strstate(ghost_task_test(accuTask)));

        ghost_task_add(accuTask);
        printf("state should be enqueued or running: %s\n",ghost_task_strstate(ghost_task_test(accuTask)));
        ghost_task_wait(accuTask);
        printf("state should be finished: %s\n",ghost_task_strstate(ghost_task_test(accuTask)));
        printf("accuTask returned %d\n",*(int *)accuTask->ret);
        ghost_task_destroy(accuTask);

        printf("\nstarting long running task @ LD0 and check if another strict LD0 task is blocking\n");
        lrTask = ghost_task_init(GHOST_TASK_FILL_LD, 0, &longRunningFunc, NULL, GHOST_TASK_DEFAULT);
        srTask = ghost_task_init(GHOST_TASK_FILL_LD, 0, &shortRunningFunc, NULL, GHOST_TASK_LD_STRICT);
        ghost_task_add(lrTask);
        ghost_task_add(srTask);

        ghost_task_waitall();
        ghost_task_destroy(lrTask);
        ghost_task_destroy(srTask);

        printf("\nstarting long running task @ LD0 and check if another non-strict LD0 task is running @ LD1\n");
        lrTask = ghost_task_init(GHOST_TASK_FILL_LD, 0, &longRunningFunc, NULL, GHOST_TASK_DEFAULT);
        srTask = ghost_task_init(GHOST_TASK_FILL_LD, 0, &shortRunningFunc, NULL, GHOST_TASK_DEFAULT);
        ghost_task_add(lrTask);
        ghost_task_add(srTask);
        
        ghost_task_waitall();

  }
	
