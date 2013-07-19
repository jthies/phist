#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

class MpiRootOnlyPrinter : public ::testing::EmptyTestEventListener 
  {
  public:
  
  // constructor
  MpiRootOnlyPrinter()
    {
    rank_=0;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
#endif
    amRoot_=(rank_==0);
    }  
     
  // Called before a test starts.
  virtual void OnTestStart(const ::testing::TestInfo& test_info) 
    {      
    if (amRoot_)
    printf("*** Test %s.%s starting.\n",             
        test_info.test_case_name(), test_info.name());    
    }
        
  // Called after a failed assertion or a SUCCEED() invocation.
  virtual void OnTestPartResult(const ::testing::TestPartResult& test_part_result) 
    {      
    if (amRoot_)
    printf("%s in %s:%d\n%s\n",
        test_part_result.failed() ? "*** Failure" : "Success",test_part_result.file_name(),
        test_part_result.line_number(),test_part_result.summary());    
    } 
  // Called after a test ends.
  virtual void OnTestEnd(const ::testing::TestInfo& test_info) 
    {
    if (amRoot_)
    printf("*** Test %s.%s ending.\n",
        test_info.test_case_name(), test_info.name());    
    }  
  
  int rank_;
  bool amRoot_;
  };
