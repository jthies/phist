#include <cpHelperFuncs.hpp>

class TYPE(CpPhistMvec) : public CpBase
{
typedef phist::ScalarTraits<_ST_> st;
private:

	TYPE(mvec_ptr) dataPtr;
	TYPE(mvec_ptr) asynData;
	int iflag;
	int update(){
		craftDbg(3, "CpPhistMvec::update \n");
		SUBR(mvec_add_mvec)(st::one(),dataPtr,st::zero(),asynData,&iflag);
		return 0;
	}

	int write( const std::string * filename){
		craftDbg(3, "CpPhistMvec: i will write now %s\n", (*filename).c_str());
		SUBR(mvec_write_bin)(asynData, (*filename).c_str() ,&iflag);  
		return 0;
	}
	
	int read(const std::string * filename){
		craftDbg(3, "CpPhistMvec: i will read now %s\n", (*filename).c_str());
		SUBR(mvec_read_bin)(asynData, (*filename).c_str() ,&iflag);   //NOTE: has to be read by asynData and then copy to dataPtr. This is because read and write operations are dependent on MPI_Comm, and dataPtr and asynData have different MPI_Comm. Thus only that vector can be read, which has written it.
		SUBR(mvec_add_mvec)(st::one(),asynData,st::zero(),dataPtr,&iflag);    
		return 0;
	}

public:
	TYPE(CpPhistMvec)(){}
	~TYPE(CpPhistMvec)(){}	

	TYPE(CpPhistMvec)( TYPE(mvec_ptr) dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD){
		dataPtr   = dataPtr_;
    cpMpiComm = cpMpiComm_;
		int nVec=0;
		iflag = 0;
		PHIST_CHK_IERR(SUBR(mvec_num_vectors)(dataPtr_,&nVec,&iflag),iflag);
		phist_const_map_ptr map;
		PHIST_CHK_IERR(SUBR(mvec_get_map)(dataPtr, &map, &iflag), iflag); 
  	PHIST_CHK_IERR(SUBR(mvec_create)(&asynData,map,nVec,&iflag),iflag);
		PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),dataPtr,st::zero(),asynData,&iflag),iflag);
	}


};

