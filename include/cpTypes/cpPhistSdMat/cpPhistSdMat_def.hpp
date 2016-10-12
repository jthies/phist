
class TYPE(cpPhistSdMat) : public cpBase
{
typedef phist::ScalarTraits<_ST_> st;
private:

	TYPE(sdMat_ptr) dataPtr;
	TYPE(sdMat_ptr) asynData;
	int iflag;

public:
	TYPE(cpPhistSdMat)(){}
	~TYPE(cpPhistSdMat)(){}	

	TYPE(cpPhistSdMat)( TYPE(sdMat_ptr) dataPtr_){
		dataPtr = dataPtr_;
		iflag = 0;
		phist_comm_ptr comm= NULL;				// TODO: comm needs to be different in the case of FT.
		PHIST_CHK_IERR(phist_comm_create(&comm, &iflag), iflag);
		int nr, nc;	
		PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(dataPtr,&nr,&iflag),iflag);
		PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(dataPtr,&nc,&iflag),iflag);
		PHIST_CHK_IERR(SUBR(sdMat_create)(&asynData, nr, nc, comm, &iflag), iflag); 

	}

	
	int update(){
		printf("cpPhistSdMat Before update: \n");
		PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),dataPtr,st::zero(),asynData,&iflag),iflag);
		return 0;
	}

	int write( const std::string * filename){
		printf("cpPhistSdMat: i will write now %s\n", (*filename).c_str());
		PHIST_CHK_IERR(SUBR(mvec_write_bin)(asynData, (*filename).c_str() ,&iflag),iflag);
		return 0;
	}
	
	int read(const std::string * filename){
		printf("cpPhistSdMat: i will read now\n");
		PHIST_CHK_IERR(SUBR(mvec_read_bin)(dataPtr, (*filename).c_str() ,&iflag),iflag);
		return 0;
	}
};
