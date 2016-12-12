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
		    MPI_Comm_set_errhandler(comm_working, errh); \
	    } \
      if ( parent != MPI_COMM_NULL || AftFailed == true)	{ \
		    if( parent != MPI_COMM_NULL){ \
			    if(*myrank==0) printf("%d_SPAWNED RANK \n", *myrank); \
			    comm_working = MPI_COMM_NULL; \
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

/*
#define AFT_RECOVER() \
  AftFailed 

#define AFT_RECOVERED() \
  AftFailed = false;

*/


/*
#define AFT_BEGIN(comm_working, myrank, argv) \
{ \
  int AftFailed = false; \
	while(success != true) \
	try{ \
     MPI_Comm parent; \
	MPI_Errhandler errh; \
	MPI_Comm_create_errhandler(&errhandler_respawn, &errh); \
	MPI_Comm_get_parent( &parent ); \
	if( parent == MPI_COMM_NULL && AftFailed == false){ \
		printf("FIRST RUNNN!!!!\n"); \
		MPI_Comm_dup (MPI_COMM_WORLD, &comm_working); \
		MPI_Comm_rank (comm_working, myrank); \
		MPI_Comm_set_errhandler(comm_working, errh); \
	} \
    if ( parent != MPI_COMM_NULL || AftFailed == true)	{ \
		if( parent != MPI_COMM_NULL){ \
			printf("%d_SPAWNED RANK \n", *myrank); \
			comm_working = MPI_COMM_NULL; \
			AftFailed = true;	\
		} \
		printf("%d: calling app_needs_repair \n", *myrank); \
		app_needs_repair(&comm_working, argv); \
		MPI_Comm_set_errhandler(comm_working, errh); \
		MPI_Comm_rank(comm_working, myrank); \
		printf("%d i am back from repair\n", *myrank); \
	} 

#define AFT_END() \
	}catch(int exception_val){ \
		sleep(1); \
		printf("exception value is %d \n", exception_val); \
		AftFailed = true;	\
 	} \
} 
*/

/*
	while(success != true)
	try{ 
     MPI_Comm parent; 
	MPI_Errhandler errh; 
	MPI_Comm_create_errhandler(&errhandler_respawn, &errh);
	MPI_Comm_get_parent( &parent ); 
	if( parent == MPI_COMM_NULL && failed == false){		// performed only in the first run 
		printf("FIRST RUNNN!!!!\n");
		MPI_Comm_dup (MPI_COMM_WORLD, &comm_working);     // for multiplezones this may need to be changed
		MPI_Comm_size (comm_working, &numprocs); 
		MPI_Comm_rank (comm_working, &myrank); 
		MPI_Comm_set_errhandler(comm_working, errh); 
	} 
    if ( parent != MPI_COMM_NULL || failed == true)	{  // spawned and remaining-alive processes. 
		if( parent != MPI_COMM_NULL){
			printf("%d/%d_SPAWNED RANK \n", myrank, numprocs); 
			comm_working = MPI_COMM_NULL; 
			failed = true;							// so that spawed also goes the restart phase. 
		}
		printf("%d/%d: calling app_needs_repair \n", myrank, numprocs); 
		app_needs_repair(&comm_working, argv, &cfg);
		MPI_Comm_rank(comm_working, &myrank);
		printf("%d i am back from repair\n", myrank); 
	}
*/



/*
	}catch(int exception_val){ 
		sleep(1); 
		printf("exception value is %d \n", exception_val);
		failed = true;	
	}
*/
