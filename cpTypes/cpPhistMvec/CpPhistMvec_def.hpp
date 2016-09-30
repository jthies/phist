
class TYPE(CpPhistMvec) : public CpBase
{
typedef phist::ScalarTraits<_ST_> st;
private:

	TYPE(mvec_ptr) dataPtr;
	TYPE(mvec_ptr) asynData;
	int iflag;

public:
	TYPE(CpPhistMvec)(){}
	~TYPE(CpPhistMvec)(){}	

	TYPE(CpPhistMvec)( TYPE(mvec_ptr) dataPtr_){
		printf("CpPhistMvec constructor is called here: \n");
		dataPtr = dataPtr_;
		int nVec;
		iflag = 0;
		PHIST_CHK_IERR(SUBR(mvec_num_vectors)(dataPtr,&nVec,&iflag),iflag);
		printf("vectors in dataPtr: %d\n",nVec);
		phist_const_map_ptr map;
		PHIST_CHK_IERR(SUBR(mvec_get_map)(dataPtr, &map, &iflag), iflag); 
  		PHIST_CHK_IERR(SUBR(mvec_create)(&asynData,map,nVec,&iflag),iflag);
		PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),dataPtr,st::zero(),asynData,&iflag),iflag);

	}

	int update(){
		printf("CpPhistMvec Before update: \n");
		PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),dataPtr,st::zero(),asynData,&iflag),iflag);
		return 0;
	}

	int write( const std::string * filename){
		printf("CpPhistMvec: i will write now %s\n", (*filename).c_str());
		PHIST_CHK_IERR(SUBR(mvec_write_bin)(asynData, (*filename).c_str() ,&iflag),iflag);
		return 0;
	}
	
	int read(const std::string * filename){
		printf("CpPhistMvec: i will read now\n");
		PHIST_CHK_IERR(SUBR(mvec_read_bin)(dataPtr, (*filename).c_str() ,&iflag),iflag);
		return 0;
	}
};

