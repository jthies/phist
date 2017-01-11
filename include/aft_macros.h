#define AFT_BEGIN(comm_working, myrank, argv) \
{ \
	int AftFailed = false; \
  do \
  { \
	  try{ \
      MPI_Comm parent; \
	    MPI_Errhandler errh; \
	    MPI_Comm_create_errhandler(&errhandler_respawn, &errh); \
	    MPI_Comm_get_parent( &parent ); \
	    if( parent == MPI_COMM_NULL && AftFailed == false){ \
		    MPI_Comm_dup (MPI_COMM_WORLD, &comm_working); \
		    MPI_Comm_rank (comm_working, myrank); \
		    if(*myrank==0) printf("FIRST RUNNN!!!!\n"); \
        getEnvParam(); \
        removeMachineFiles(&comm_working); \
        initRescueNodeList(&comm_working); \
		    MPI_Comm_set_errhandler(comm_working, errh); \
	    } \
      if ( parent != MPI_COMM_NULL || AftFailed == true)	{ \
		    if( parent != MPI_COMM_NULL){ \
			    if(*myrank==0) printf("%d_SPAWNED RANK \n", *myrank); \
			    comm_working = MPI_COMM_NULL; \
          getEnvParam(); \
			    AftFailed = true;	\
		    } \
		    if(*myrank==0) printf("%d: calling app_needs_repair \n", *myrank); \
		    app_needs_repair(&comm_working, argv); \
		    MPI_Comm_set_errhandler(comm_working, errh); \
		    MPI_Comm_rank(comm_working, myrank); \
		    if(*myrank==0) printf("%d i am back from repair\n", *myrank); \
      }

#define AFT_END() \
    AftFailed = false; \
	  }catch(int exception_val){ \
		  sleep(1); \
		  printf("exception value is %d \n", exception_val); \
		  AftFailed = true;	\
 	  } \
  } \
  while(AftFailed == true);\
} 


