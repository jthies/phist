/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
  TEST_F(CLASSNAME, get_comm) 
  {
    phist_const_comm_ptr comm;
    phist_map_get_comm(map_,&comm,&iflag_);
    ASSERT_EQ(0,iflag_);
  }


  /*! Test the comm_get_rank function. */
  TEST_F(CLASSNAME, maps_compatible_with_same_map)
  {
    if (problemTooSmall_) return;
    phist_maps_compatible(map_,map_,&iflag_);
    ASSERT_EQ(0,iflag_);
    if (defaultMap_!=map_)
    {
      phist_maps_compatible(map_,map_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    // create the default map again and see if it is compatible with the orginal one
    phist_map_ptr tmp_map=NULL;
    phist_map_create(&tmp_map,comm_,_N_,&iflag_);
    ASSERT_EQ(iflag_,0);
    phist_maps_compatible(defaultMap_,tmp_map,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_map_delete(tmp_map,&iflag_);
    ASSERT_EQ(0,iflag_);
  }


  /*! Test the comm_get_rank function. */
  TEST_F(CLASSNAME, defaultMap_is_linear) 
  {
    if (problemTooSmall_) return;
    phist_const_comm_ptr comm;
    phist_map_get_comm(map_,&comm,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    int rank, nproc;
    phist_comm_get_rank(comm,&rank,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_comm_get_size(comm,&nproc,&iflag_);
    ASSERT_EQ(0,iflag_);
  
    phist_gidx ilower, iupper;
    phist_map_get_ilower(defaultMap_,&ilower,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_map_get_iupper(defaultMap_,&iupper,&iflag_);
    ASSERT_EQ(0,iflag_);
    // make sure the map is monotone
    ASSERT_TRUE(ilower<=iupper);
    // make sure the map is 0-based and the last GID is _N_-1
    ASSERT_TRUE(rank!=0 || ilower==0);
    ASSERT_TRUE(rank!=nproc-1 || iupper==_N_-1);
    
    // make sure the local length of the map is iupper-ilower+1
    phist_lidx nloc;
    phist_map_get_local_length(defaultMap_,&nloc,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    ASSERT_EQ(nloc,iupper-ilower+1);
    
#ifdef PHIST_HAVE_MPI
    // check if sum(nloc)=nglob
    int i_nloc, i_nglob;
    i_nloc=(int)nloc;
    iflag_=MPI_Allreduce(&i_nloc,&i_nglob,1,MPI_INT,MPI_SUM,mpi_comm_);
    EXPECT_EQ((int)MPI_SUCCESS,iflag_);
    ASSERT_EQ(i_nglob,(int)_N_);
    // gather all the offsets and check they are monotonous and contiguous
    int i_ilower=(int)ilower,i_iupper=(int)iupper;
    int i_all_ilower[nproc], i_all_iupper[nproc];
    iflag_=MPI_Allgather(&i_ilower,1,MPI_INT,i_all_ilower,1,MPI_INT,mpi_comm_);
    EXPECT_EQ((int)MPI_SUCCESS,iflag_);
    bool contig=true;
    for (int i=0; i<nproc-1; i++)
    {
      contig&=(i_all_iupper[i]=i_all_ilower[i+1]);
    }
#endif  
  }

