
class TYPE(CpPhistSdMat) : public CpBase
{
typedef phist::ScalarTraits<_ST_> st;
private:

	TYPE(sdMat_ptr) dataPtr;
	TYPE(sdMat_ptr) asynData;
	int iflag;

public:
	TYPE(CpPhistSdMat)(){}
	~TYPE(CpPhistSdMat)(){}	

	TYPE(CpPhistSdMat)( TYPE(sdMat_ptr) dataPtr_, const MPI_Comm cpMpiComm_=MPI_COMM_WORLD ){
		dataPtr = dataPtr_;
		iflag = 0;
    cpMpiComm = cpMpiComm_;
		phist_comm_ptr comm= NULL;				// TODO: comm needs to be different in the case of FT.
		PHIST_CHK_IERR(phist_comm_create(&comm, &iflag), iflag);
		int nr, nc;	
		PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(dataPtr,&nr,&iflag),iflag);
		PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(dataPtr,&nc,&iflag),iflag);
		PHIST_CHK_IERR(SUBR(sdMat_create)(&asynData, nr, nc, comm, &iflag), iflag); 
	}

	
	int update(){
		craftDbg(3, "CpPhistSdMat Before update: \n");
    SUBR(sdMat_add_sdMat)(st::one(),dataPtr,st::zero(),asynData,&iflag);
		return 0;
	}

	int write( const std::string * filename){
		craftDbg(3, "CpPhistSdMat: i will write now %s\n", (*filename).c_str());
		SUBR(mvec_write_bin)(asynData, (*filename).c_str() ,&iflag);
		return 0;
	}
	
	int read(const std::string * filename){
		craftDbg(3, "CpPhistSdMat: i will read now\n");
		SUBR(mvec_read_bin)(asynData, (*filename).c_str() ,&iflag);   //NOTE: has to be read by asynData and then copy to dataPtr. This is because read and write operations are dependent on MPI_Comm, and dataPtr and asynData have different MPI_Comm. Thus only that vector can be read, which has written it.
    SUBR(sdMat_add_sdMat)(st::one(),asynData,st::zero(),dataPtr,&iflag);
		return 0;
	}
};
