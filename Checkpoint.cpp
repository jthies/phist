#include "helperFuncs.hpp"
#include "Checkpoint.hpp"
#include "CpPOD.hpp"

void Checkpoint::update(){
	for(cp_const_map:: iterator it = objects.begin(); it != objects.end(); ++it)
	{
		std::cout << "updating: " << it->first << std::endl;
		it->second->update();
	}	


}

void Checkpoint::write()
{
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = NumberToString(myrank_);
	std::string filename;
  	for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
  	{
		filename = cpPath + "/" + (it->first) + myrank + ".ckpt"; 
    	std::cout << it->first << "Checkpoint::write is called here " << filename <<std::endl;
    	it->second->write(&filename);
  	}
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
    	std::cout << it->first << "Checkpoint::read is called here "<<std::endl;
    	it->second->read(&filename);
  	}
}


