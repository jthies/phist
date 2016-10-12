#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <openssl/md5.h>

#include "include/checkpoint.hpp"



Checkpoint::Checkpoint(const std::string cpBasePath_=exec("pwd"), const MPI_Comm cpMpiComm_=MPI_COMM_WORLD){
	static int idx=0;										
	++idx;	
//	name 				= name_;
	name 				= SSTR(idx);
	cpMpiComm 	= cpMpiComm_;	
	cpBasePath 	= cpBasePath_;
	cpPath 			= cpBasePath_ + "/" + "Checkpoint-" + name;
//	std::cout << "cpPath:" << cpPath << std::endl;
	cpIdx				= idx;
	cpVersionPrefix = "v-";
	cpCommitted = false;
	cpVersion 	= 0;
	numBufferCps=2;
//	std::cout << "constructor is being called here idx = " << name << std::endl;
	if(cpMpiComm == MPI_COMM_WORLD) 
		std::cout << "MPI_COMM_WORLD it is \n" << std::endl;
	mkCpDir(cpPath);	
#ifdef SCR																		// check if CPAFTLIB was compiled with SCR
	useSCR 		= true;
#else
	useSCR 		= false;
#endif

}	

Checkpoint::~Checkpoint(){
#ifdef SCR
	if(useSCR == true){
		printf("====== SCR_Finalize is done ====== \n");
		SCR_Finalize();
	}
#endif
}

void Checkpoint::disableSCR(){
#ifdef SCR
	printf("SCR is disabled. Use SCR to reduce the impact of Checkpoint/restart\n");
	useSCR = false;
#endif
}

void Checkpoint::commit(){
	cpCommitted = true;	
}

int Checkpoint::mkCpDir(std::string path){
	std::string cmd = "mkdir -p " + path;
	system (cmd.c_str());
	return EXIT_SUCCESS;
}

int Checkpoint::writeCpMetaData(){			// writes a tmp file first and then replaces it after complete write
	int myrank = -1;	
	MPI_Comm_rank(cpMpiComm, &myrank);
	if(myrank == 0){
		std::string filename, filenameTmp;
		filenameTmp = cpPath + "/" + "metadata.ckpt.tmp"; 
		filename = cpPath + "/" + "metadata.ckpt"; 
		std::ofstream fstr;
		fstr.open ((filenameTmp).c_str(), std::ios::out );	
		if(fstr.is_open()){
			fstr << cpVersion << std::endl;
			fstr.close();
		}
		std::string cmd = "mv " + filenameTmp + " " + filename;
		system (cmd.c_str());
	}
	return EXIT_SUCCESS;
}

int Checkpoint::readCpMetaData(){
	std::string filename, line;
	filename = cpPath + "/" + "metadata.ckpt"; 
	std::ifstream fstr;
	fstr.open ((filename).c_str(), std::ios::in );	
	if(fstr.is_open()){
		while ( getline (fstr,line) )
		{
			cpVersion = atoi(line.c_str());
		}
		fstr.close();
		std::cout << "cpVersion: " << cpVersion << std::endl;
	}
	return EXIT_SUCCESS;
}
int Checkpoint::deleteBackupCp(){
	int myrank = -1;	
	MPI_Comm_rank(cpMpiComm, &myrank);
	if(myrank == 0){
		std::string toRmDir = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion-numBufferCps-1);
		std::string cmd = "rm -r " + toRmDir;	
		struct stat sb;
		if (stat(toRmDir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
		{
			system ( cmd.c_str());
 		}
	}
	return EXIT_SUCCESS;
}


int Checkpoint::update(){
	for(cp_const_map:: iterator it = objects.begin(); it != objects.end(); ++it)
	{
		std::cout << "updating: " << it->first << std::endl;
		it->second->update();
	}	
	return EXIT_SUCCESS;
}

int Checkpoint::write()				// TODO: make two version of write. 1) PFS 2) SCR. and do MPI/IO for PFS 
{					
#ifdef SCR
		SCRwrite();
#else 
		PFSwrite();
#endif
	return EXIT_SUCCESS;
}

int Checkpoint::SCRwrite(){
#ifdef SCR
	if(useSCR == true){
		SCR_Start_checkpoint();
	}
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion);

  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
		char * fNameScrPath = new char[256]; 
		if(useSCR == true){
			char * tmpFilename = new char[256];
			strcpy(tmpFilename, filename.c_str()); 
			SCR_Route_file(tmpFilename, fNameScrPath);
			printf("new filename %s\n", fNameScrPath);
			filename = fNameScrPath;
		}
//    	std::cout << it->first << " Checkpoint:write is called here " << filename <<std::endl;
    it->second->write(&filename);
	}
	if(useSCR == true){
		int valid = 1;
		SCR_Complete_checkpoint(valid);     // valid flag should be 1 if every CP is successfully written. TODO: check the flag from each write call.
	}
#endif
	return EXIT_SUCCESS;
}


int Checkpoint::PFSwrite(){
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion);
	mkCpDir(cpPathVersion);	
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{
#ifdef MPIIO																										// Parallel PFS IO
		filename = cpPathVersion + "/" + (it->first) + ".ckpt"; 
#else																														// Serial PFS IO
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
#endif
		if(EXIT_SUCCESS != it->second->write(&filename))	
				return EXIT_FAILURE;
	}
	MPI_Barrier(cpMpiComm);			// TODO: do MPI_Gather here 
	writeCpMetaData();
	++cpVersion;
	deleteBackupCp();
	return EXIT_SUCCESS;
}


int Checkpoint::read()				// TODO: make two version of read. 1) PFS 2) SCR
{
#ifdef SCR
		SCRread();
#else 
		PFSread();
#endif
	return EXIT_SUCCESS;
}

int Checkpoint::PFSread(){
	readCpMetaData();									// cpVersion is updated here
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion); 
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{			
#ifdef MPIIO																										// Parallel PFS IO 
		filename = cpPathVersion + "/" + (it->first) + ".ckpt";				
#else																														// Serial PFS IO 
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
#endif
//    	std::cout << it->first << " Checkpoint::read is called here "<<std::endl;
		if(EXIT_SUCCESS != it->second->read(&filename))	
			return EXIT_FAILURE;
  }
	MPI_Barrier(cpMpiComm);
	++cpVersion;
	return EXIT_SUCCESS;
}

int Checkpoint::SCRread(){
#ifdef SCR
	readCpMetaData();						// cpVersion is updated here
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion); 
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{			
		filename = cpPathVersion + "/" + (it->first) + ".ckpt"; 
		char * fNameScrPath = new char[256]; 
		char * tmpFilename = new char[256];
		strcpy(tmpFilename, filename.c_str()); 
		SCR_Route_file(tmpFilename, fNameScrPath);
		printf("new filename %s\n", fNameScrPath);
		filename = fNameScrPath;
//    	std::cout << it->first << " Checkpoint::read is called here "<<std::endl;
		if(EXIT_SUCCESS != it->second->read(&filename))	
			return EXIT_FAILURE;
  }
	MPI_Barrier(cpMpiComm);
	++cpVersion;
#endif

	return EXIT_SUCCESS;
}

