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
        getEnvParam(); \
        removeMachineFiles(&comm_working); \
        initRescueNodeList(&comm_working); \
		    MPI_Comm_set_errhandler(comm_working, errh); \
	    } \
      if ( parent != MPI_COMM_NULL || AftFailed == true)	{ \
		    if( parent != MPI_COMM_NULL){ \
			    comm_working = MPI_COMM_NULL; \
          getEnvParam(); \
			    AftFailed = true;	\
		    } \
		    craftDbg(1, "%d: calling app_needs_repair ", *myrank); \
		    app_needs_repair(&comm_working, argv); \
		    MPI_Comm_set_errhandler(comm_working, errh); \
		    MPI_Comm_rank(comm_working, myrank); \
        craftTime("repairTime", &comm_working);\
		    craftDbg(1, "%d i am back from repair", *myrank); \
      }

#define AFT_END() \
    AftFailed = false; \
	  }catch(int exception_val){ \
      craftTime("failTime");\
		  AftFailed = true;	\
 	  } \
  } \
  while(AftFailed == true);\
} 


