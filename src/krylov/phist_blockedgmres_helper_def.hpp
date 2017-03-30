/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
// small class for a ring buffer for the subspaces
// TODO: use a real template instead of this "TYPE"-stuff?
class TYPE(MvecRingBuffer)
{
  public:
    TYPE(MvecRingBuffer)(int size) : mvecs_(size,NULL), mvecs_used_(size,0), lastIndex_(0) {}

    // we can handle failures probably more cleanly if not done in the constructor
    void create_mvecs(phist_const_map_ptr map, int nvecs, int* iflag)
    {
      for(int i = 0; i < mvecs_.size(); i++)
      {
        PHIST_CHK_IERR(*iflag = (mvecs_[i] != NULL) ? -1 : 0, *iflag);
        PHIST_CHK_IERR(SUBR( mvec_create ) (&mvecs_[i], map, nvecs, iflag), *iflag);
      }
    }

    // must be called before the destructor
    void delete_mvecs(int *iflag)
    {
      for(int i = 0; i < mvecs_.size(); i++)
      {
        PHIST_CHK_IERR(SUBR( mvec_delete ) (mvecs_[i], iflag), *iflag);
        mvecs_[i] = NULL;
      }
    }

    ~TYPE(MvecRingBuffer)()
    {
      int iflag = 0;
      // print an error if there are still allocated mvecs
      // as this could use up all memory quite fast if used multiple times!
      for(int i = 0; i < mvecs_.size(); i++)
      {
        PHIST_CHK_IERR(iflag = (mvecs_[i] != NULL) ? -1 : 0, iflag);
      }
    }

    // get vector at index i
    TYPE(mvec_ptr) at(int i) {return mvecs_.at(i);}

    int size() {return mvecs_.size();}

    // just to make sure no element is used twice
    void incRef(int i) {mvecs_used_.at(i)++;}
    void decRef(int i) {mvecs_used_.at(i)--;}
    int refCount(int index) const {return mvecs_used_.at(index);}

    // get last used index
    int lastIndex() {return lastIndex_;}

    // get preceeding index
    int prevIndex(int i, int n = 1) {return (i+size()-n)%size();}

    // get next index and make sure it is unused
    void getNextUnused(int &nextIndex, int *iflag)
    {
      *iflag = 0;
      nextIndex = (lastIndex_+1)%mvecs_.size();
      if( mvecs_used_[nextIndex] != 0 )
      {
        *iflag = -1;
        return;
      }
      lastIndex_ = nextIndex;
    }

  private:
    std::vector<TYPE(mvec_ptr)> mvecs_;
    std::vector<int> mvecs_used_;
    int lastIndex_;

    // hide copy constructor etc
    TYPE(MvecRingBuffer)(const TYPE(MvecRingBuffer)&);
    const TYPE(MvecRingBuffer)& operator=(const TYPE(MvecRingBuffer)&);
};




