#include "Checkpoint.hpp"
#include <openssl/md5.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int Checkpoint::update(){
	for(cp_const_map:: iterator it = objects.begin(); it != objects.end(); ++it)
	{
		std::cout << "updating: " << it->first << std::endl;
		it->second->update();
	}	
	return EXIT_SUCCESS;
}

int Checkpoint::write()				// TODO: make two version of write. 1) PFS 2) SCR
{					
#ifdef SCR
	if(useSCR == true){
		SCR_Start_checkpoint();
	}
#endif
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion);
#ifndef SCR
	mkCpDir(cpPathVersion);	
#endif

  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
#ifdef SCR
		char * fNameScrPath = new char[256]; 
		if(useSCR == true){
			char * tmpFilename = new char[256];
			strcpy(tmpFilename, filename.c_str()); 
			SCR_Route_file(tmpFilename, fNameScrPath);
			printf("new filename %s\n", fNameScrPath);
			filename = fNameScrPath;
		}
#endif	
//    	std::cout << it->first << " Checkpoint:write is called here " << filename <<std::endl;
    it->second->write(&filename);
	}
	
#ifndef SCR
	MPI_Barrier(cpMpiComm);
	writeCpMetaData();
	++cpVersion;
	deleteBackupCp();
	MPI_Barrier(cpMpiComm);
#endif	
	
#ifdef SCR
	if(useSCR == true){
		int valid = 1;
		SCR_Complete_checkpoint(valid);     // valid flag should be 1 if every CP is successfully written. TODO: check the flag from each write call.
	}
#endif
	return EXIT_SUCCESS;
}

int Checkpoint::read()				// TODO: make two version of read. 1) PFS 2) SCR
{
	readCpMetaData();
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion); 
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{			
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
#ifdef SCR 
		char * fNameScrPath = new char[256]; 
		if(useSCR == true){
			char * tmpFilename = new char[256];
			strcpy(tmpFilename, filename.c_str()); 
			SCR_Route_file(tmpFilename, fNameScrPath);
			printf("new filename %s\n", fNameScrPath);
			filename = fNameScrPath;
		}
#endif
//    	std::cout << it->first << " Checkpoint::read is called here "<<std::endl;
		if(EXIT_SUCCESS != it->second->read(&filename))	
			return EXIT_FAILURE;
  }
	MPI_Barrier(cpMpiComm);
	++cpVersion;
	return EXIT_SUCCESS;
}

