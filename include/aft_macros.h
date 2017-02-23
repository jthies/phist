#define AFT_BEGIN(CRAFT_comm_working, CRAFT_myrank, CRAFT_argv) \
{ \
	int CRAFT_aftFailed = false; \
  do \
  { \
	  try{ \
      MPI_Comm CRAFT_parent; \
	    MPI_Errhandler CRAFT_errh; \
	    MPI_Comm_create_errhandler(&CRAFT_errhandlerRespawn, &CRAFT_errh); \
	    MPI_Comm_get_parent( &CRAFT_parent ); \
	    if( CRAFT_parent == MPI_COMM_NULL && CRAFT_aftFailed == false){ \
		    MPI_Comm_dup (MPI_COMM_WORLD, &CRAFT_comm_working); \
		    MPI_Comm_rank (CRAFT_comm_working, CRAFT_myrank); \
        CRAFT_getEnvParam(); \
        CRAFT_removeMachineFiles(&CRAFT_comm_working); \
        CRAFT_initRescueNodeList(&CRAFT_comm_working); \
		    MPI_Comm_set_errhandler(CRAFT_comm_working, CRAFT_errh); \
	    } \
      if ( CRAFT_parent != MPI_COMM_NULL || CRAFT_aftFailed == true)	{ \
		    if( CRAFT_parent != MPI_COMM_NULL){ \
			    CRAFT_comm_working = MPI_COMM_NULL; \
          CRAFT_getEnvParam(); \
			    CRAFT_aftFailed = true;	\
		    } \
		    craftDbg(1, "%d: calling CRAFT_appNeedsRepair ", *CRAFT_myrank); \
		    CRAFT_appNeedsRepair(&CRAFT_comm_working, CRAFT_argv); \
		    MPI_Comm_set_errhandler(CRAFT_comm_working, CRAFT_errh); \
		    MPI_Comm_rank(CRAFT_comm_working, CRAFT_myrank); \
        craftTime("repairTime", &CRAFT_comm_working);\
		    craftDbg(1, "%d i am back from repair", *CRAFT_myrank); \
      }

#define AFT_END() \
    CRAFT_aftFailed = false; \
	  }catch(int CRAFT_exception_val){ \
      craftTime("failTime");\
		  CRAFT_aftFailed = true;	\
 	  } \
  } \
  while(CRAFT_aftFailed == true);\
} 


