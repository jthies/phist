#include "include/aft.h"
#include <algorithm>

static char estr[MPI_MAX_ERROR_STRING]=""; 
static int strl; /* error messages */

int initRescueNodeList(MPI_Comm * const comm){    // executed once for the FIRST RUN.
  int myrank, numprocs;
  MPI_Comm_rank(*comm, &myrank);
  MPI_Comm_size(*comm, &numprocs);
  std::string line = "";
  MPI_Barrier(*comm);
  std::string filenameP = pbsNodeFile;
  std::string filenameF = machinefileFailedProcs;
  std::string filenameA = machinefileActiveProcs;
  
  std::ifstream fstrA, fstrP;
  std::set<std::string> pbsNodeList;
  // ====== PBS_NODEFILE List ====== // 
  if(myrank==getFirstRank(comm)){
	  fstrP.open ((filenameP).c_str(), std::ios::in );	
	  if(fstrP.is_open()){
      while (getline (fstrP, line))
      {
        pbsNodeList.insert(line); 
      }
		  fstrP.close();
    }else{
		  std::cerr << "Can't open file " << filenameP << std::endl;			
		  return EXIT_FAILURE;
	  }
  }
 
  // ====== CREATE ACTIVE NODE LIST ===== // 
  std::vector<std::string> activeMachineList;
  writeActiveMachineList(&activeMachineList, comm);
  MPI_Barrier(*comm);
  std::set<std::string> activeNodeList;
  copyVecToSet(&activeMachineList, &activeNodeList);

  // ====== RESCUE = PBSLIST - ACTIVELIST ====== // This is only true at the start of the job. Afterwards, the Rescue list has to be managed.

  std::set<std::string> rescueNodeList;
  findSetDiffence(&pbsNodeList, &activeNodeList, &rescueNodeList);
  if(myrank==getFirstRank(comm)){
    printSet("PBS_NODELIST", pbsNodeList);  
    printSet("ACTIVE_NODELIST", activeNodeList);  
    printSet("RESCUE_NODELIST", rescueNodeList);  
  }
  writeSetList(machinefileRescueProcs, rescueNodeList, comm);

  MPI_Barrier(*comm);
  sync();
  return EXIT_SUCCESS;
}

int writeActiveMachineList(std::vector<std::string> * activeMachineList_, MPI_Comm * const comm){
  makeActiveMachineList(activeMachineList_, comm);
  writeVectorList(machinefileActiveProcs, *activeMachineList_, comm);
  return 0;
}

int makeActiveMachineList(std::vector<std::string> *activeMachineList_, MPI_Comm * const comm){
  int numprocs, myrank;
	MPI_Comm_size(*comm, &numprocs);
	MPI_Comm_rank(*comm, &myrank);
	
	char *node_name= (char *) malloc(sizeof(char)*MPI_MAX_PROCESSOR_NAME);
	int len;
	MPI_Get_processor_name(node_name, &len);
  char ** activeMachineList=NULL;
		activeMachineList 		= (char ** ) malloc (sizeof(char *)*numprocs);
		for(int i=0; i<numprocs ; ++i){
			activeMachineList[i] 			= (char * ) malloc (sizeof(char)*len);
		}
		
	for(int i=0;i<numprocs; ++i){
		sprintf(activeMachineList[i] , "%s", "");
	}
	if(myrank!=0){
		MPI_Send(node_name, len+1, MPI_CHAR, 0, 0, *comm);
	}
	
	MPI_Status * my_recv_status = (MPI_Status *) malloc (sizeof(MPI_Status)*numprocs);
	char *node_name_rcvd= (char *) malloc(sizeof(char)*len);
	sprintf(node_name_rcvd, "%s", "");
	if(myrank==getFirstRank(comm)){
		sprintf(activeMachineList[0], "%s", node_name);
		for(int i=1; i<numprocs; ++i){
			sprintf(node_name_rcvd, "%s", "");
			MPI_Recv(node_name_rcvd, len+1, MPI_CHAR, i, 0, *comm, &my_recv_status[i]);
			sprintf(activeMachineList[i], "%s", node_name_rcvd);
		}
	}

	MPI_Barrier(*comm);

  for(int i=0;i<numprocs; ++i){
    activeMachineList_->push_back(activeMachineList[i]);
  }

	MPI_Barrier(*comm);
  return 0;
}


/* Do all the magic in the error handler i.e. MPIX_Comm_revoke */
void errhandler_respawn(MPI_Comm* pcomm, int* errcode, ...) {

    int eclass;
    MPI_Error_class(*errcode, &eclass);
    int myrank = -1;
    MPI_Comm_rank(*pcomm, &myrank);
        
    MPI_Error_string(*errcode, estr, &strl);
    craftDbg(1, "%04d: errhandler invoked with error %s\n", myrank, estr);

//    fprintf(stderr, "%04d: errhandler invoked with error %s\n", myrank, estr);
	
    if( MPIX_ERR_PROC_FAILED != eclass &&
        MPIX_ERR_REVOKED != eclass ) {
        MPI_Abort(MPI_COMM_WORLD, *errcode);
    }
    MPIX_Comm_revoke(*pcomm);
    craftDbg(1, "Comm is removked. Now throwing error.");
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
	craftDbg(1, "appNeedsRepair start");
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
//          app_reload_ckpt(world);
//    	}

  if(craftCommRecoveryPolicy == "SHRINKING") {
  	craftDbg(3, "Initiating SHRINKING recovery.");
    craftTime("before shrink");
		MPIX_Comm_shrink(*comm, tempcomm);
    craftTime("after shrink");
  	craftDbg(3, "SHRINKING recovery done.");
  }
  if(craftCommRecoveryPolicy == "NON-SHRINKING") {
  	craftDbg(3, "Initiating NON-SHRINKING recovery.");
  	MPIX_Comm_replace(*comm, tempcomm , argv);
  	craftDbg(3, "NON-SHRINKING recovery done.");
  }
	MPI_Barrier(*tempcomm);
  *comm = *tempcomm;
//	MPI_Comm_dup(*tempcomm, comm);
	MPI_Comm_rank(*comm, &myrank);	
	craftDbg(1, "appNeedsRepair done");
  return 0; /* we have repaired the world, we need to reexecute */
}

int MPIX_Comm_replace(MPI_Comm comm, MPI_Comm *newcomm, char** argv) {
  MPI_Comm icomm, // the intercomm between the spawnees and the old (shrinked) world 
           scomm, // the local comm for each sides of icomm : shrinked comm
   	       mcomm; // the intracomm, merged from icomm 
  MPI_Group cgrp, sgrp, dgrp;
  int rc, flag, rflag, i, nc, ns, nd, crank=-1, srank=-1, drank=-1;
	 
redo:
  if (comm == MPI_COMM_NULL){ // am I a newly spawned process?
    // I am a new spawnee, waiting for my new rank assignment
    // it will be sent by rank 0 in the old world 
		MPI_Comm_get_parent(&icomm);
		scomm = MPI_COMM_WORLD;
		MPI_Recv(&crank, 1, MPI_INT, 0, 1, icomm, MPI_STATUS_IGNORE);
  }
	else {
    // kill_all_procs_on_failed_processhost(comm , my_param);		// TODO: fix this thing. after killing the process, no other process detects the error,... make it fault tolerant as well.
    // I am a survivor: Spawn the appropriate number
    // of replacement processes (we check that this operation worked
    // before we procees further) 
    // First: remove dead processes 
    craftTime("before shrink");
		MPIX_Comm_shrink(comm, &scomm);
    craftTime("after shrink");
    
	  MPI_Comm_rank(scomm, &srank );
	  MPI_Comm_set_errhandler(scomm, MPI_ERRORS_RETURN);
    MPI_Barrier(scomm);
    nd = getNumFailed(&comm, &scomm);
    int * failedRanks= new int[nd];
    for( int i = 0; i < nd; ++i){ failedRanks[i] = -1; }
    getFailedRanks(&comm, &scomm, failedRanks);    // updates nd and failedRanks 
    for(int i=0; i < nd; ++i){
      craftDbg(1, "failed_rank[%d] is %d", i, failedRanks[i]);
//      craftLog(comm, "===== This is the %d entry of LOG  =====\n", tempx);
	  }
	  
   // ====== TODO: Determine the spawnHosts list based on 'spawn policy' here =====
	 //We handle failures during this fuction execution ourseves//
    std::string spawnHostsTmp; 
    craftTime("before makeSpawnList");
    makeSpawnList(nd, failedRanks, &scomm, &spawnHostsTmp);     // This function makes list where to spawn, depending on the spawn policy (REUSE, NO-REUSE)
    craftTime("after makeSpawnList");
	  char * scr_copy_location = new char[256];
	  sprintf(scr_copy_location, "PARTNER");
    // ======= Spawn Processes ====== // 
    char * machinefileSpawnTmp = new char[1024];
    sprintf(machinefileSpawnTmp, "%s", machinefileSpawnProcs);
	  MPI_Info spawn_info;
	  MPI_Info_create(&spawn_info);
 
//    MPI_Info_set(spawn_info, "hostfile", machinefileSpawnTmp);  // writing host file creates problems sometimes due to buffer IO.
    char * spawnHosts = new char[spawnHostsTmp.size()];
    sprintf(spawnHosts, "%s", spawnHostsTmp.c_str());
    MPI_Info_set(spawn_info, "host", spawnHosts);
    MPI_Info_set(spawn_info, "SCR_COPY_TYPE", scr_copy_location);
    craftTime("before spawning");
	  rc = MPI_Comm_spawn(argv[0], argv+1, nd, spawn_info, 0, scomm, &icomm, MPI_ERRCODES_IGNORE);
    craftTime("after spawning");
	  flag = (MPI_SUCCESS == rc);
    craftTime("before agree");
	  MPIX_Comm_agree(scomm, &flag);
    craftTime("after agree");
	  if (!flag) {	// spawed has failed
		  if ( MPI_SUCCESS == rc ) {
			  MPIX_Comm_revoke(icomm);
			  MPI_Comm_free(&icomm);
		  }
		  MPI_Comm_free(&scomm);
		  craftDbg(1, "%04d: comm_spawn failed, redo", crank);
		  goto redo;
	  }
	 
   // ======= Spawn Processes ====== // 
 
//	  remembering the former rank: we will reassign the same ranks in the new world 
    craftTime("before MPI_Send");
    MPI_Comm_rank (comm, &crank);
	  MPI_Comm_rank (scomm, &srank);
	  // the rank 0 in the scomm comm is going to determine the ranks at which the spares need to be inserted. 
	  // determining the ranks of the dead processes 
	  if (0 == srank) {
		  for(int i=0; i<nd; i++) {
        drank = failedRanks[i];
			  craftDbg(2, "Rank of dead process is (drank): %d", drank);
			  MPI_Send(&drank, 1, MPI_INT, i, 1, icomm);
		  }
	  }	
    craftTime("after MPI_Send");
	  }
 //   Merge the intercomm, to reconstruct an intracomm ( we check 
 //	 that this operation worked before we proceeed further 
	  rc = MPI_Intercomm_merge(icomm, 1 , &mcomm);
    craftTime("after merge", &mcomm);
	  rflag = flag = (MPI_SUCCESS == rc);
    craftTime("before AGREE 1", &mcomm);
	  MPIX_Comm_agree(scomm, &flag);
    craftTime("after AGREE 1", &mcomm);
	  if( MPI_COMM_WORLD != scomm) MPI_Comm_free(&scomm);
    craftTime("before AGREE 2", &mcomm);
	  MPIX_Comm_agree(icomm, &rflag);
    craftTime("after AGREE 2", &mcomm);
	  MPI_Comm_free(&icomm);
	  if(!(flag && rflag) ){
		  if(MPI_SUCCESS == rc){
			  MPI_Comm_free(&mcomm);
		  }
		  craftDbg(1, "%04d: Intercomm_merge failed, redo", crank);
		  goto redo;
	  }
    craftTime("after 2 agrees", &mcomm);
	  int myrank_mcomm= -1;
	  MPI_Comm_rank(mcomm, &myrank_mcomm);
	  MPI_Barrier(mcomm);
	  craftDbg(2, "%d_Merge is done", myrank_mcomm);
 //   Merge is done. Now, reorder the mcomm according to the original rank ordering in comm
 //	 * Split does the magic. removing spare processes and reordering ranks so that all 
 //	 * surviving processes remain at their former places
    craftTime("before split", &mcomm);
 	  rc = MPI_Comm_split(mcomm, 1, crank, newcomm);
    craftTime("after split", newcomm );
 //   Split or some of the communications above may have failed if 
//	 new failures have disrupted the process, we need to make sure 
 //	 we succeeded at all the ranks, or retry until it works. 
	  flag = (MPI_SUCCESS==rc);
	  MPIX_Comm_agree(mcomm, &flag);
	  MPI_Comm_free(&mcomm);
	  if( !flag ) {
		  if (MPI_SUCCESS == rc) {
			  MPI_Comm_free (newcomm);
		  }
		  craftDbg(1, "%04d: comm_split failed, redo", crank);
		  goto redo;
	  }	
    std::vector<std::string> activeNodeList;
    craftTime("before write active machine list", newcomm );
    writeActiveMachineList(&activeNodeList, newcomm); 
    craftTime("after write active machine list", newcomm );
	  MPI_Barrier(*newcomm);
	  craftDbg(3, "%d_Split is done", myrank_mcomm);
  
//	 restore the error handler 
	  if (MPI_COMM_NULL != comm){
		  MPI_Errhandler errh;
		  MPI_Comm_get_errhandler (comm, &errh);
		  MPI_Comm_set_errhandler (*newcomm, errh);
	  }
    craftTime("replace function finished", newcomm);
//    printNodeName(newcomm);
    return MPI_SUCCESS;
 }


int printNodeName( const MPI_Comm * const comm){
	 int numprocs, myrank;
	 MPI_Comm_size(*comm, &numprocs);
	 MPI_Comm_rank(*comm, &myrank);
	 
	 char *node_name= (char *) malloc(sizeof(char)*MPI_MAX_PROCESSOR_NAME);
	 int len;
	 MPI_Get_processor_name(node_name, &len);
	 for(int i=0;i<numprocs; ++i){
		 if(myrank==i){
       printf("====== myrank:%d, node_name:%s\n", myrank, node_name);	
     }
     MPI_Barrier(*comm);
   }
  return 0;
}

int getNumFailed(const MPI_Comm * const comm, const MPI_Comm * scomm){
  int nd = 0;
  MPI_Group cGrp, sGrp, dGrp;
	MPI_Comm_group(*comm, &cGrp);
	MPI_Comm_group(*scomm, &sGrp);
	MPI_Group_difference(cGrp, sGrp, &dGrp);
  MPI_Group_size(dGrp, &nd);
  if( nd == 0 ) {
		  craftDbg(1, "No process failure was detected in MPIX_Comm_shrink ( nd=%d )", nd);
			return -1;
	}
  return nd;
}

int getFailedRanks(const MPI_Comm * const comm, const MPI_Comm * scomm, int * failedRanks){
  // comm = original broken comm, scomm = shrinked comm , dGrp = dead group.
  MPI_Group cGrp, sGrp, dGrp;
	MPI_Comm_group(*comm, &cGrp);
	MPI_Comm_group(*scomm, &sGrp);
	MPI_Group_difference(cGrp, sGrp, &dGrp);
  int nd;
  MPI_Group_size(dGrp, &nd);
	for(int i=0; i < nd; ++i){
		MPI_Group_translate_ranks(dGrp, 1, &i, cGrp, &failedRanks[i]);
 		craftDbg(1, "failed_rank[%d] is %d", i, failedRanks[i]);
	}

  return EXIT_SUCCESS;
}

int makeSpawnList(const int nd, const int * const failedRanks, MPI_Comm * const comm, std::string * spawnList){
  MPI_Barrier (*comm);
  if(craftCommSpawnPolicy == "REUSE"){
    makeSpawnListReuse(nd, failedRanks, comm, spawnList);
  }
  if(craftCommSpawnPolicy == "NO-REUSE"){
//    writeFailedList(nd, failedRanks, comm);         // This writing is not necessary for functionality but just for keeping record of the failed nodes.
    makeSpawnListNoReuse(nd, failedRanks, comm, spawnList);   // also updates/writes the list of Rescue Processes
  }
  MPI_Barrier(*comm); 
  sync();
//  if(craftCommSpawnPolicy == "DYNAMIC"){
//    craftDbg(3, "craftCommSpawnPolicy == DYNAMIC");   
//  }
  return 0;
}

int makeSpawnListReuse(const int nd, const int * const failedRanks, MPI_Comm * const comm, std::string * spawnList){
  craftDbg(3, "craftCommSpawnPolicy == REUSE");   
  int myrank=-1;
  MPI_Comm_rank(*comm, &myrank);
  std::string tempStr;
  std::string line;
  std::ifstream fstri;
  std::string filenamei = machinefileActiveProcs;

	fstri.open ((filenamei).c_str(), std::ios::in );	
	if( fstri.is_open() ){
    int i=0;
    while ( getline (fstri,line) )
    {
      for(int j=0;j<nd;++j){
        if (failedRanks[j] == i){   // failedRanks matches an entry in machinefileActiveProcs 
          tempStr += line + ",";
        }
      }
      ++i;
    }
		fstri.close();
    tempStr.erase(tempStr.end()-1);
    *spawnList = tempStr;
    craftDbg(1, "spawnList reuse %s", spawnList->c_str());
//    if(myrank == getFirstRank(comm)) std::cout << "spawnList reuse= " << *spawnList << '\n';
  }
  else{
		if(myrank == getFirstRank(comm)) std::cerr << "Can't open file " << filenamei << std::endl;			
		return EXIT_FAILURE;
	}
  return 0;
}

int makeSpawnListNoReuse(const int nd, const int * const failedRanks, MPI_Comm * const comm, std::string * spawnList){
  // Opens the machinefileRescueProcs 
  // and writes the first node name in spawnList string.
  // renew the file with remaining available nodes.
  // TODO: Due to bug in MPI_Comm_spawn, processes can only be spawned on one node. 
  // Thus, only one first entry of rescue-nodelist file need to be written in spawn file. 

  craftDbg(3, "makeSpawnListNoReuse ");
  int myrank=-1;
  MPI_Comm_rank(*comm, &myrank);
//  ===== READ FULL RESCUE FILE ===== // 
  std::set<std::string> rescueNodeList;
  std::string line;
  std::string filenameR = machinefileRescueProcs;
  std::ifstream fstrR;
	fstrR.open ((filenameR).c_str(), std::ios::in );	
	if(fstrR.is_open() ){
    while ( getline (fstrR,line) )
    { 
      rescueNodeList.insert(line);
    }
  }else{
		  if(myrank == getFirstRank(comm)) std::cerr << "Can't open file " << filenameR << std::endl;			
		  return EXIT_FAILURE;
	}
  if(rescueNodeList.empty()){
    makeSpawnListReuse(nd, failedRanks, comm, spawnList); // if no rescue node is available //  
  } 
  else{
    int nd_ = nd;   // TODO: remove if above limitation is removed. 
    nd_ = 1;
    std::set<std::string>::iterator it;
    it = rescueNodeList.begin();
    for(int i = 0; i < nd_; ++i)
    { 
      *spawnList = *it + ",";
      rescueNodeList.erase(it);      // updating the rescue Node list for future.
    }
    printSet("== NEW RESCUE SET NO-REUSE ==", rescueNodeList); 
    writeSetList(machinefileRescueProcs, rescueNodeList, comm);
  }
    
  MPI_Barrier(*comm);
  return 0;
}

int writeFailedList(const int nd, const int * const failedRanks, MPI_Comm * const comm){
  craftDbg(3, "craftCommSpawnPolicy == NO-REUSE");   
  int myrank=-1;
  MPI_Comm_rank(*comm, &myrank);
  if(myrank == getFirstRank(comm))
  {
    std::string line;
    std::ifstream fstri;
    std::ofstream fstro;
    std::string filenamei = machinefileActiveProcs;
    std::string filenameo = machinefileFailedProcs;
    std::set<std::string> failedNodeNames;
  
	  fstri.open ((filenamei).c_str(), std::ios::in );	
	  if(fstri.is_open()){
      int i=0;
      while ( getline (fstri,line) )
      {
        for(int j=0;j<nd;++j){
          if (failedRanks[j] == i){   // failedRanks matches an entry in machinefileActiveProcs 
            failedNodeNames.insert(line);
          }
        }
        ++i;
      }
		  fstri.close();
    }
    else{
		  std::cerr << "Can't open file " << filenamei << "\n or can't open file "<< filenameo << std::endl;			
		  return EXIT_FAILURE;
	  }

	  fstro.open ((filenameo).c_str(), std::ios::app );	
    if(fstro.is_open()){
      std::set<std::string>::iterator it;
      for ( it=failedNodeNames.begin(); it != failedNodeNames.end(); ++it )
      {
            fstro <<  *it << std::endl;
      }
		  fstro.close();
    }
  }
  MPI_Barrier(*comm);
  return 0;
}

int removeMachineFiles(MPI_Comm * const comm){
  removeFile( machinefileActiveProcs, comm);
  removeFile( machinefileFailedProcs, comm);
  removeFile( machinefileSpawnProcs , comm);
  removeFile( craftLogFile , comm);
  return 0;
}

int removeFile(const char * filename, MPI_Comm * const comm){
  int myrank = -1;	
	MPI_Comm_rank(*comm, &myrank);
	if(myrank == getFirstRank(comm)){
    char * cmd = new char[1024];
    sprintf( cmd , "rm %s", filename);
		system ( cmd );
	} 
  return 0;
}

int getFirstRank(MPI_Comm* const comm){
  int rank, size;
  MPI_Comm_rank(*comm, &rank);
  MPI_Comm_size(*comm, &size);
  int * ranks = new int[size];
  MPI_Allgather(&rank, 1, MPI_INT, ranks, 1, MPI_INT, *comm);
  return ranks[0];
}

int findSetDiffence(std::set<std::string> * s1, std::set<std::string> * s2, std::set<std::string> * res){
  std::set_difference(s1->begin(), s1->end(), s2->begin(), s2->end(), std::inserter(*res, res->begin())); 
  return 0;
}

int printSet(std::string toPrint, std::set<std::string> s){
  std::set<std::string>::iterator it;
  craftDbg(4, "==== %s ====", toPrint.c_str());
  for ( it= s.begin(); it != s.end(); ++it){
    craftDbg(4, "%s", (*it).c_str());
  }
  return 0;
}

int writeSetList(const std::string setFileName, std::set<std::string> set, MPI_Comm * const comm){
  int myrank(-1);
  MPI_Comm_rank(*comm, &myrank); 
  if(myrank==getFirstRank(comm)){
    std::set<std::string>::iterator it;
    std::ofstream fstr;  
    printSet(setFileName.c_str(), set); 

    // ===== WRITE RESCUE LIST ===== // 
    fstr.open ((setFileName).c_str(), std::ios::out );	
    if(fstr.is_open()){
      for(it=set.begin(); it != set.end() ; ++it){
        fstr << *it << std::endl;
      }
      fstr.close(); 
      sync();
    }else{
		  std::cerr << "Can't open file " << setFileName << std::endl;			
		  return EXIT_FAILURE;
	  }
  }
  return EXIT_SUCCESS;
}

int writeVectorList(const std::string vecFileName, std::vector<std::string> vec, MPI_Comm * const comm){
  int myrank(-1);
  MPI_Comm_rank(*comm, &myrank); 
  if(myrank==getFirstRank(comm)){
    std::vector<std::string>::iterator it;
    std::ofstream fstr;  
//    printVector(vecFileName.c_str(), vec); 
    // ===== WRITE RESCUE LIST ===== // 
    fstr.open ((vecFileName).c_str(), std::ios::out );	
    if(fstr.is_open()){
      for(it=vec.begin(); it != vec.end() ; ++it){
        fstr << *it << std::endl;
      }
      fstr.close(); 
      sync();
    }else{
		  std::cerr << "Can't open file " << vecFileName << std::endl;			
		  return EXIT_FAILURE;
	  }
  }
  return EXIT_SUCCESS;
}


int copyVecToSet(std::vector<std::string> * vecSrc, std::set<std::string> *setDst){
  std::vector<std::string>::iterator it;
  for (it = vecSrc->begin(); it != vecSrc->end() ; ++it){
    setDst->insert(*it);    
  }
  return EXIT_SUCCESS;
}



int craftTimeToFile(const std::string fileName, MPI_Comm * const comm){
  double t=0.0;
  int myrank=-1;
  MPI_Comm_rank(*comm, &myrank);
  get_walltime_ (&t);
  if(myrank == getFirstRank(comm))
  {
    std::ofstream fstr;
    fstr.open ((fileName).c_str(), std::ios::out );	
    if(fstr.is_open()){
      fstr.precision(16);
      fstr << t << std::endl;
      fstr.close(); 
      sync();
    }else{
		  std::cerr << "Can't open file " << fileName << std::endl;			
		  return EXIT_FAILURE;
	  } 
  }
  return EXIT_SUCCESS;
}

int craftTimeToFile(const std::string fileName){
 double t=0.0;
 int myrank=-1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  get_walltime_ (&t);
  if(myrank == 0)
  {
    std::ofstream fstr;
    fstr.open ((fileName).c_str(), std::ios::out );	
    if(fstr.is_open()){
      fstr.precision(16);
      fstr << t << std::endl;
      fstr.close(); 
      sync();
    }else{
		  std::cerr << "Can't open file " << fileName << std::endl;			
		  return EXIT_FAILURE;
	  } 
  }
  return EXIT_SUCCESS;
}


