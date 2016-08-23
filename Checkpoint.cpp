#include "Checkpoint.hpp"

void Checkpoint::update(){
	for(cp_const_map:: iterator it = objects.begin(); it != objects.end(); ++it)
	{
		std::cout << "updating: " << it->first << std::endl;
		it->second->update();
	}	


}

void Checkpoint::write()
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
  	for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
  	{
			filename = cpPath + "/" + (it->first) + myrank + ".ckpt"; 
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
    	std::cout << it->first << " Checkpoint:write is called here " << filename <<std::endl;
    	it->second->write(&filename);
  	
		}

#ifdef SCR
	if(useSCR == true){
		int valid = 1;
		SCR_Complete_checkpoint(valid);     // valid flag should be 1 if every CP is successfully written. TODO: check the flag from each write call.
	}
#endif

}

void Checkpoint::read()
{
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
  	for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
  	{
			filename = cpPath + "/" + (it->first) + myrank + ".ckpt"; 
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
    	std::cout << it->first << "Checkpoint::read is called here "<<std::endl;
    	it->second->read(&filename);
  	}
}




				


