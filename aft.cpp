#include "aft.h"

static char estr[MPI_MAX_ERROR_STRING]=""; static int strl; /* error messages */


/* Do all the magic in the error handler */
void errhandler_respawn(MPI_Comm* pcomm, int* errcode, ...) {

    int eclass;
    MPI_Error_class(*errcode, &eclass);
    int myrank = -1;
    MPI_Comm_rank(*pcomm, &myrank);
        
    MPI_Error_string(*errcode, estr, &strl);
    fprintf(stderr, "%04d: errhandler invoked with error %s\n", myrank, estr);
	
    if( MPIX_ERR_PROC_FAILED != eclass &&
        MPIX_ERR_REVOKED != eclass ) {
        MPI_Abort(MPI_COMM_WORLD, *errcode);
    }
    MPIX_Comm_revoke(*pcomm);
    printf("%d_Before throwing the error \n", myrank);
    throw 5;
}

/* repair comm world, reload checkpoints, etc...
 *  Return: true: the app needs to redo some iterations
 *  false: no failure was fixed, we do not need to redo any work.
 */
int app_needs_repair(MPI_Comm *comm, char ** argv) {
    /* This is the first time we see an error on this comm, do the swap of the
 *      * worlds. Next time we will have nothing to do. */
	MPI_Comm *tempcomm = new MPI_Comm[1];
	int myrank=-1;
	if(*comm != MPI_COMM_NULL){
		MPI_Comm_rank(*comm, &myrank);
	}
//	if( *comm != MPI_COMM_NULL ) {
	        /* swap the worlds */
        /* We keep comm around so that the error handler remains attached until the
         * user has completed all pending ops; it is expected that the user will
         * complete all ops on comm before posting new ops in the new world.
         * Beware that if the user does not complete all ops on comm and the handler
         * is invoked on the new world inbetween, comm may be freed while
         * operations are still pending on it, and a fatal error may be
         * triggered when these ops are finally completed (possibly in Finalize)*/
//        	if( MPI_COMM_NULL != world ) MPI_Comm_free(&world);
//        	if( world == comm ) return false; /* ok, we repaired nothing, no need to redo any work */
	 //     app_reload_ckpt(world);
//    	}

	printf("%d_calling MPIX_Comm_replace now \n", myrank);
	MPIX_Comm_replace(*comm, tempcomm , argv);
	printf("%d_MPIX_Comm_replace is done \n", myrank);
	
	MPI_Barrier(*tempcomm);
//	*comm = *tempcomm;
	MPI_Comm_dup(*tempcomm, comm);
	MPI_Comm_rank(*comm, &myrank);	
	printf("%d_=========== MPIX_Comm_replace is done  ============ \n", myrank);
    return 0; /* we have repaired the world, we need to reexecute */
}



int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm, char** argv) {
	static int failCount = 0;
    	int verbose = true;
	char * machinefile = new char[256];
	sprintf (machinefile , "machinefile");
	MPI_Comm icomm, /* the intercomm between the spawnees and the old (shrinked) world */
	            scomm, /* the local comm for each sides of icomm : shrinked comm*/
     	       mcomm; /* the intracomm, merged from icomm */
     MPI_Group cgrp, sgrp, dgrp;
     int rc, flag, rflag, i, nc, ns, nd, crank=-1, srank=-1, drank=-1;
     static int num_procs_failed=0;
	static int total_num_procs_failed=0;
	
redo:
	if (comm == MPI_COMM_NULL){ // am I a newly spawned process?
		/* I am a new spawnee, waiting for my new rank assignment
 		* it will be sent by rank 0 in the old world */
		printf("I am the spawned process \n");
		MPI_Comm_get_parent(&icomm);
		scomm = MPI_COMM_WORLD;
		MPI_Recv(&crank, 1, MPI_INT, 0, 1, icomm, MPI_STATUS_IGNORE);
	}
	else {
// 			kill_all_procs_on_failed_processhost(comm , my_param);		// TODO: fix this thing. after killing the process, no other process detects the error,... make it fault tolerant as well.
 

		/* I am a survivor: Spawn the appropriate number
		* of replacement processes (we check that this operation worked
		* before we procees further) */
		/* First: remove dead processes */
		MPIX_Comm_shrink(comm, &scomm);
		MPI_Comm_size(scomm, &ns);
		MPI_Comm_size(comm, &nc);
		nd = nc - ns; 				/*number of dead processes*/	
		total_num_procs_failed += nd; 
		if( nd == 0 ) {
			if( verbose ) {
				printf("No process failure was detected in MPIX_Comm_shrink ( nd=%d )\n", nd);
			}
			MPI_Comm_free(&scomm);
			*newcomm = comm;
			return MPI_SUCCESS;
		}
	/*We handle failures during this fuction execution ourseves*/
	MPI_Comm_set_errhandler(scomm, MPI_ERRORS_RETURN);
	MPI_Comm_rank(scomm, &srank );
//	if(srank==0){
//		generate_host_file_spawn( nc, nd, total_num_procs_failed, my_param);
//		printf("......generating host file is done\n");
//	}
			
	char * spawnHosts = new char[256]; 
	sprintf(spawnHosts, "");
	makeHostList( machinefile, failCount, nd, nc, spawnHosts);
	printf("calc. spawnHosts is %s \n", spawnHosts);

	char * scr_copy_location = new char[256];
	sprintf(scr_copy_location, "PARTNER");
	MPI_Info spawn_info;
	MPI_Info_create(&spawn_info);
//	MPI_Info_set(spawn_info, "hostfile", "host_file_spawn");
	MPI_Info_set(spawn_info, "host", spawnHosts);
	MPI_Info_set(spawn_info, "SCR_COPY_TYPE", scr_copy_location);

	rc = MPI_Comm_spawn(argv[0], argv+1, nd, spawn_info, 0, scomm, &icomm, MPI_ERRCODES_IGNORE);
	flag = (MPI_SUCCESS == rc);
	MPIX_Comm_agree(scomm, &flag);
	if (!flag) {	/*spawed has failed*/
		if ( MPI_SUCCESS == rc ) {
			MPIX_Comm_revoke(icomm);
			MPI_Comm_free(&icomm);
		}
		MPI_Comm_free(&scomm);
		if (verbose) fprintf(stderr, "%04d: comm_spawn failed, redo\n", crank);
		goto redo;
	}
	failCount += nd;	
	
	/* remembering the former rank: we will reassign the same ranks in the new world */
	MPI_Comm_rank (comm, &crank);
	MPI_Comm_rank (scomm, &srank);
	/* the rank 0 in the scomm comm is going to determine the ranks at which the spares need to be inserted. */
	/* determining the ranks of the dead processes */
	if (0 == srank) {
		/*getting the gruop of dead processes
		 * those in comm, but not in scomm are the dead processes */
		MPI_Comm_group(comm, &cgrp);
		MPI_Comm_group(scomm, &sgrp);
		MPI_Group_difference(cgrp, sgrp, &dgrp);
		/* Computing the rank assignment for the newly inserted spares */
		for(i=0; i<nd; i++) {
			MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
			if (verbose){
				printf("Rank of dead process is (drank): %d\n", drank);
			}
			/* sending the calculated rank to the spawned process */
			MPI_Send(&drank, 1, MPI_INT, i, 1, icomm);
		}
		MPI_Group_free(&cgrp); 
		MPI_Group_free(&sgrp); 
		MPI_Group_free(&dgrp); 
	}	
	}
	/* Merge the intercomm, to reconstruct an intracomm ( we check 
	that this operation worked before we proceeed further */	
	
	rc = MPI_Intercomm_merge(icomm, 1 , &mcomm);
	rflag = flag = (MPI_SUCCESS == rc);
	MPIX_Comm_agree(scomm, &flag);
	if( MPI_COMM_WORLD != scomm) MPI_Comm_free(&scomm);
	MPIX_Comm_agree(icomm, &rflag);
	MPI_Comm_free(&icomm);
	if(!(flag && rflag) ){
		if(MPI_SUCCESS == rc){
			MPI_Comm_free(&mcomm);
		}
		if(verbose); fprintf(stderr, "%04d: Intercomm_merge failed, redo\n", crank);
		goto redo;
	}
	int myrank_mcomm= -1;
	MPI_Comm_rank(mcomm, &myrank_mcomm);
	MPI_Barrier(mcomm);
	if(myrank_mcomm == 0) 	printf("%d_Merge is done\n", myrank_mcomm);
	/* Merge is done. Now, reorder the mcomm according to the original rank ordering in comm
	* Split does the magic. removing spare processes and reordering ranks so that all 
	* surviving processes remain at their former places*/
 	rc = MPI_Comm_split(mcomm, 1, crank, newcomm);
	/* Split or some of the communications above may have failed if 
	new failures have disrupted the process, we need to make sure 
	we succeeded at all the ranks, or retry until it works. */
	flag = (MPI_SUCCESS==rc);
	MPIX_Comm_agree(mcomm, &flag);
	MPI_Comm_free(&mcomm);
	if( !flag ) {
		if (MPI_SUCCESS == rc) {
			MPI_Comm_free (newcomm);
		}
		if (verbose) fprintf(stderr, "%04d: comm_split failed, redo\n", crank);
		goto redo;
	}	
	MPI_Barrier(*newcomm);
	if(myrank_mcomm == 0)	printf("%d_Split is done\n", myrank_mcomm);
	sleep(1);

	/* restore the error handler */	
	if (MPI_COMM_NULL != comm){
		MPI_Errhandler errh;
		MPI_Comm_get_errhandler (comm, &errh);
		MPI_Comm_set_errhandler (*newcomm, errh);
	}
    	printf("%d Failures so far %d \n", crank, total_num_procs_failed);
//    	cerr<<"("<<crank<<")Failures so far :"<< total_num_procs_failed <<endl;
     return MPI_SUCCESS;
}

void makeHostList(const char * machinefile, const int failCount, const int nd, const int nc, char * spawnHosts){
	int counter = nc + failCount;
	FILE * fp_read;
	if(NULL == (fp_read =fopen(machinefile,"r")) ){
		fprintf(stderr, "Error: Unable to open file (%s)\n",machinefile);
	}
	int i=0, j=0;
	char * tempVar = new char[256];
	sprintf(tempVar, "");
	for(j=0;j<(nc + failCount);++j){
		fscanf(fp_read, "%s", tempVar);
	}
	for(i=(nc + failCount) ; i < (nc + failCount + nd ) ; ++i){
		fscanf(fp_read, "%s", tempVar);
		if((nc + failCount + nd)-1 > i) // a comma should separate between host names
		{
			sprintf(tempVar, "%s,", tempVar);
//			printf("tempVar %s\n", tempVar);
		}
		sprintf(spawnHosts, "%s%s", spawnHosts, tempVar);
//		printf("spawnHosts %s\n", spawnHosts);
	}
	fclose(fp_read);
	printf("spawnHosts %s\n", spawnHosts);
	return ;
}

/*
int generate_host_file_spawn(int nc,  int num_procs_failed, int total_num_procs_failed, ConfigFile * cfg){

//	====== determine which nodes to start the process (host) ===== only done by master process =====	 //

	char * temp_var 		= (char *) malloc (sizeof(char)*256);
	char * host_file_spawn 	= (char *) malloc (sizeof(char)*256);
	sprintf(host_file_spawn, "%s", "host_file_spawn");
	FILE * fp_read;
	if(NULL == (fp_read =fopen(m->machine_file,"r")) ){
		fprintf(stderr, "Error: Unable to open file (%s)\n",cfg->machinefile);
	}
	FILE * fp_write;
	if(NULL == (fp_write =fopen(host_file_spawn,"w+")) ){
		fprintf(stderr, "Error: Unable to open file (%s)\n",host_file_spawn);
	}
	int i=0, j=0;
	for(j=0;j<(nc + total_num_procs_failed - num_procs_failed);++j){
		fscanf(fp_read, "%s", temp_var);
	}
	for(i=(nc + total_num_procs_failed - num_procs_failed) ; i < (nc + total_num_procs_failed) ; ++i){
		fscanf(fp_read, "%s", temp_var);
		fprintf(fp_write, "%s\n", temp_var);
		printf("new spawned host %s\n", temp_var);
	}
	fclose(fp_read);
	fclose(fp_write);
	return 0;
}
*/

/*
int make_active_machine_list(MPI_Comm * comm_working){ // create a global list of active nodes. First every node sends its name to 0. 0th-process then gives this list to everyone else.
	int numprocs, myrank, i, j; 
	MPI_Comm_size(*comm_working, &numprocs);
	MPI_Comm_rank(*comm_working, &myrank);
	
	char *node_name= (char *) malloc(sizeof(char)*MPI_MAX_PROCESSOR_NAME);
	int len;
	MPI_Get_processor_name(node_name, &len);
	
	if(active_machine_list==NULL)
	{
		active_machine_list 		= (char ** ) malloc (sizeof(char *)*numprocs);
		for(i=0; i<numprocs ; ++i){
			active_machine_list[i] 			= (char * ) malloc (sizeof(char)*len);
		}
	}
		
	for(i=0;i<numprocs; ++i){
		sprintf(my_param->active_machine_list[i] , "%s", "");
	}
	if(myrank!=0){
		MPI_Send(node_name, len+1, MPI_CHAR, 0, 0, *comm_working);
	}
	
	MPI_Status * my_recv_status = (MPI_Status *) malloc (sizeof(MPI_Status)*numprocs);
	char *node_name_rcvd= (char *) malloc(sizeof(char)*len);
	sprintf(node_name_rcvd, "%s", "");
	if(myrank==0){
		sprintf(my_param->active_machine_list[0], "%s", node_name);
		for(i=1; i<numprocs; ++i){
			sprintf(node_name_rcvd, "%s", "");
			MPI_Recv(node_name_rcvd, len+1, MPI_CHAR, i, 0, *comm_working, &my_recv_status[i]);
			sprintf(my_param->active_machine_list[i], "%s", node_name_rcvd);
		}
	}
	MPI_Barrier(*comm_working);
	// give the list to every other process now
	// NOTE: maybe it is much faster to write a file than to do MPI Send Recv in this manner
	if(myrank==0){	// send to others
		char *node_name_send= (char *) malloc(sizeof(char)*len);
		for(j=1;j<numprocs;++j){
			for(i=0;i<numprocs;++i){
				sprintf(node_name_send, "%s", my_param->active_machine_list[i]);
				MPI_Send(node_name_send, len+1, MPI_CHAR, j, 0, *comm_working);

			}
		}
	}
	
	if(myrank!=0){	// recv the node-names
		for(i=0;i<numprocs;++i){
			MPI_Recv(node_name_rcvd, len+1, MPI_CHAR, 0, 0, *comm_working, &my_recv_status[i]);
			sprintf(my_param->active_machine_list[i], "%s", node_name_rcvd);
		}
	}
	for(i=0;i<numprocs; ++i){
		if(myrank==i){
			for(j=0; j<numprocs; ++j){
				printf("%d: node_name_rcvd %s\n", myrank, my_param->active_machine_list[j]);
			}
		}
	}
	return 0;
}
*/



/*
int kill_all_procs_on_failed_processhost(MPI_Comm * comm , param * my_param ){
	
/// get ranks of failed process, get their machine(host) names from my_param->active_machine_list. compare myrank with my_param->active_machine_list. if they match kill the process. 

	MPI_Group failed_grp;
	MPI_Group comm_working_grp;
	MPI_Comm_group(*comm, &comm_working_grp);
	int myrank;
	MPI_Comm_rank(*comm, &myrank);
// 	MPIX_Comm_agree(comm, &flag);				// Agree before failure_ack?
	MPIX_Comm_failure_ack(*comm);
	MPIX_Comm_failure_get_acked(*comm, &failed_grp);			// All processes are aware of the failures 

	int failed_grp_size;
	MPI_Group_size(failed_grp, &failed_grp_size);
// 	printf("failed_grp_size is %d\n", failed_grp_size);
	int i,j, failed_rank=0;
	int * failed_procs_list = (int *) malloc (sizeof(int)*failed_grp_size);
	for(i=0; i< failed_grp_size; ++i){
		MPI_Group_translate_ranks(failed_grp, 1, &i, comm_working_grp, &failed_procs_list[i]);
// 		printf("failed_rank[%d] is %d\n", i, failed_procs_list[i]);
	}
	
	char ** failed_procs_machine_list = (char **) malloc (sizeof(char *)*failed_grp_size);
	for(i=0;i<failed_grp_size;++i){
		failed_procs_machine_list[i] = (char *) malloc (sizeof(char)*MPI_MAX_PROCESSOR_NAME);
	}
	
	for(i=0; i<failed_grp_size; ++i){
		sprintf(failed_procs_machine_list[i],"%s",  my_param->active_machine_list[ failed_procs_list[i] ]);
// 		printf("failed_procs_machine_list[%d]: %s\n", i, failed_procs_machine_list[i]);
	}
	
	
	char *my_node_name= (char * )malloc (sizeof(char)*MPI_MAX_PROCESSOR_NAME);
	int len;
	MPI_Get_processor_name(my_node_name, &len);
// 	printf("my_node_name: %s\n", my_node_name);
	for(i=0; i<failed_grp_size; ++i){
		if (strcmp(failed_procs_machine_list[i],my_node_name) == 0)
		{
			printf("%d: ... SIGKILL NOW\n", myrank);
			kill(getpid(), SIGKILL);
		}
		if (strcmp(failed_procs_machine_list[i],my_node_name) != 0)
		{
			printf("%d: ....... not killing,  because failed_procs_machine_list[i] %s, my_node_name %s\n", myrank, failed_procs_machine_list[i], my_node_name);
		}
	}
	MPI_Barrier(*comm);
	MPIX_Comm_revoke(*comm);
// 	int flag; 
// 	MPIX_Comm_agree(*comm, );
	return 0;
}
*/

/*
int MPIX_Comm_replace(MPI_Comm worldwspares, MPI_Comm comm, MPI_Comm *newcomm) {
    MPI_Comm shrinked;
    MPI_Group cgrp, sgrp, dgrp;
    int rc, flag, i, nc, ns, nd, crank, srank, drank;
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("%d: comm_replace is called \n", myrank);
    while(1){
    	// First: remove dead processes 
    	OMPI_Comm_shrink(worldwspares, &shrinked);
		printf("shrink done\n");
	// We do not want to crash if new failures come... 
    	MPI_Comm_set_errhandler( shrinked, MPI_ERRORS_RETURN );
    	MPI_Comm_size(shrinked, &ns); MPI_Comm_rank(shrinked, &srank);

		MPI_Comm_group(shrinked, &sgrp);
		int sgrp_size=0, sgrp_rank=0;
		MPI_Group_size(sgrp, &sgrp_size);
		MPI_Group_rank(sgrp, &sgrp_rank);
		printf("rank %d of sgrp of size %d\n", sgrp_rank, sgrp_size);
		
		
    	if(MPI_COMM_NULL != comm) { // I was not a spare before... 
        	// not enough processes to continue, aborting. 
        	MPI_Comm_size(comm, &nc); if( nc > ns ) MPI_Abort(comm, MPI_ERR_INTERN);
        	// remembering the former rank: we will reassign the same
          	 // ranks in the new world. 
        	MPI_Comm_rank(comm, &crank);

        	// the rank 0 in the shrinked comm is going to determine the 
          	 // ranks at which the spares need to be inserted. 
       	 	if(0 == srank) {
            		// getting the group of dead processes: 
                 	//   those in comm, but not in shrinked are the deads 
            		MPI_Comm_group(comm, &cgrp); MPI_Comm_group(shrinked, &sgrp);
	            	MPI_Group_difference(cgrp, sgrp, &dgrp); MPI_Group_size(dgrp, &nd);
					printf("========== Number of dead processes %d \n", nd);
        	    	// Computing the rank assignment for the newly inserted spares 
        	    	for(i=0; i<ns-(nc-nd); i++) {
               		 	if( i < nd ) MPI_Group_translate_ranks(dgrp, 1, &i, cgrp, &drank);
                		else drank=-1; // still a spare //
	                	// sending their new assignment to all spares 
        	        	MPI_Send(&drank, 1, MPI_INT, i+nc-nd, 1, shrinked);
						printf("drank is: %d\n", drank);
            		}
            		MPI_Group_free(&cgrp); MPI_Group_free(&sgrp); MPI_Group_free(&dgrp);
        	}
    	} 
	else { // I was a spare, waiting for my new assignment //
        	MPI_Recv(&crank, 1, MPI_INT, 0, 1, shrinked, MPI_STATUS_IGNORE);
    	}

    	// Split does the magic: removing spare processes and reordering ranks
         /// so that all surviving processes remain at their former place 
    	rc = MPI_Comm_split(shrinked, crank<0?MPI_UNDEFINED:1, crank, newcomm);

    	// Split or some of the communications above may have failed if  
	 // new failures have disrupted the process: we need to 
 	 // make sure we succeeded at all ranks, or retry until it works. 
    	flag = (MPI_SUCCESS==rc);
    	OMPI_Comm_agree(shrinked, &flag);
    	MPI_Comm_free(&shrinked);
    	if( !flag ) {
        	MPI_Comm_free( newcomm );
    	}
    	else
        	break;
    }
    return MPI_SUCCESS;
}
*/



