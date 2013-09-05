#include <iostream>
#include "ghost.h"
#include "ghost_util.h"
#include "gtest/gtest.h"

//! This class suppresses output from any MPI process but rank 0 in googletest,
//! both to screen and an XML file. In order to activate it, create an instance
//! and append it to the gtest listeners. It then takes control of the default-
//! listeners and removes them from the list.
class MpiRootOnlyPrinter : public ::testing::EmptyTestEventListener 
  {
  public:
  
  // constructor
  MpiRootOnlyPrinter()
    {
    rank_=ghost_getRank(MPI_COMM_WORLD);
    amRoot_=(rank_==0);
    // take ownership of the default printer
    myPrinter_= 
    ::testing::UnitTest::GetInstance()->listeners().Release(
    ::testing::UnitTest::GetInstance()->listeners().default_result_printer()
    );
    myXmlPrinter_= 
    ::testing::UnitTest::GetInstance()->listeners().Release(
    ::testing::UnitTest::GetInstance()->listeners().default_xml_generator()
    );
    // become the default printer
    }
    
  ~MpiRootOnlyPrinter()
    {
    delete myPrinter_;
    if (myXmlPrinter_!=NULL) delete myXmlPrinter_;
    }

  virtual void OnTestCaseStart(const ::testing::TestCase& test_case) 
    {
    if (amRoot_) myPrinter_->OnTestCaseStart(test_case);
    }

  // Called before a test starts.
  virtual void OnTestStart(const ::testing::TestInfo& test_info) 
    {
    if (amRoot_) myPrinter_->OnTestStart(test_info);
    }
        
  // Called after a failed assertion or a SUCCEED() invocation.
  virtual void OnTestPartResult(const ::testing::TestPartResult& test_part_result) 
    {
    if (amRoot_) myPrinter_->OnTestPartResult(test_part_result);
    } 
  // Called after a test ends.
  virtual void OnTestEnd(const ::testing::TestInfo& test_info) 
    {
    if (amRoot_) myPrinter_->OnTestEnd(test_info);
    }  
    
void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration)
  {
    if (amRoot_) 
      {
      myPrinter_->OnTestIterationEnd(unit_test,iteration);
      if (myXmlPrinter_!=NULL) myXmlPrinter_->OnTestIterationEnd(unit_test,iteration);
  }   }

  
  ::testing::TestEventListener* myPrinter_;
  ::testing::TestEventListener* myXmlPrinter_;
  int rank_;
  bool amRoot_;
  };
