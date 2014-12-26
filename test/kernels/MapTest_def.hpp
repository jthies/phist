#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithMap<_N_> {
public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    iflag_=0;
    KernelTestWithMap<_N_>::SetUp();
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithMap<_N_>::TearDown();
    }

};

  /*! Test the comm_get_rank function. */
  TEST_F(CLASSNAME, get_comm) {
        const_comm_ptr_t comm;
        phist_map_get_comm(map_,&comm,&iflag_);
	ASSERT_EQ(0,iflag_);
}

