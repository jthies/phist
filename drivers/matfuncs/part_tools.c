#include "matfuncs.h"

ghost_gidx_t perm2d( ghost_gidx_t row, gidx_t arg2 )
{
  static int n1=-1, n2;
  static int np, np1, np2, pid, pid1, pid2;
  static int n1_loc, n2_loc;
  static int off1, off2;
  static ghost_gidx_t global_offset;
  
  // first call: initialize
  if (n1==-1)
  {
    n1 = row;
    n2 = arg2;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
  
    np1=1;
    np2=np;
    while (np1<=np)
    {
      n1_loc = n1/np1;
      n2_loc = n2/np2;
      if (n1_loc==n2_loc) break;
      np1++;
      np2=np/np1;
    }
    if (np1*np2!=np || n1_loc!=n2_loc || np1*n1_loc!=n1 || np2*n2_loc!=n2)
    {
        PHIST_SOUT(PHIST_WARNING,"our simple cartesian partitioning is not\n"
                                 "sophisticated enough for this setup, reverting\n"
                                 "to standard linear partitioing\n");
        n1=-2;
    }
    else
    {
      pid1 = pid%np1;
      pid2 = (pid-pid1)/np1;
      global_offset = pid*(n1_loc*n2_loc);
      off1=pid1*n1_loc;
      off2=pid2*n2_loc;
    }
    PHIST_SOUT(PHIST_VERBOSE,"2D partitioning: %d x %d domains of size %d x %d\n",
          np1, np2, n1_loc, n2_loc);
  }
  else if (n1==-2)
  {
    // partitioning failed
    return row;
  }
  else if (arg2==-1)
  {
    // iperm: convert new GID to original one. This is used
    // to convert the given row index, so it is always a
    // local variable.
    int lrow=(int)(row-global_offset);
    int i1 = lrow%n1_loc;
    int i2 = (lrow-i1)/n1_loc;
    return (ghost_gidx_t)(off2 + i2)*n1+off1+i1;
  }
  else if (arg2==+1)
  {
    // perm: convert old GID to new one

    int i1=row%n1;      // global indices (i1,i2)
    int i2=(row-i1)/n1;

    // to which process does it belong now?
    int p1 = i1/n1_loc;
    int p2 = i2/n2_loc;
    int p  = p2*np1+p1;
    ghost_gidx_t offset = p*(ghost_gidx_t)(n1_loc*n2_loc);

    // local indices (j1,j2) in new partitioning
    int j1=i1%n1_loc;
    int j2=i2%n2_loc;
    return offset+j2*n1_loc+j1;
  }
  return -1; // should never get here  
}


ghost_gidx_t perm3d( ghost_gidx_t row, gidx_t arg2, gidx_t arg3 )
{
  static int n1=-1, n2, n3;
  static int np, np1, np2, np3, pid, pid1, pid2, pid3;
  static int n1_loc, n2_loc, n3_loc;
  static int off1, off2, off3;
  static ghost_gidx_t global_offset;
  
  // first call: initialize
  if (n1==-1)
  {
    n1 = row;
    n2 = arg2;
    n3 = arg3;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    np1=1;
    np2=1;
    np3=np;
    // we assume that the grid can be partitioned
    // into np cube-shaped subdomains, so for instance
    // if np=6 we allow (n1,n2,n3)=(100,200,300) but not (100,100,100).
    // The latter works, however, with np=8.
    while (np1<=np && np2<=np && np3!=0)
    {
      n1_loc = n1/np1;
      n2_loc = n2/np2;
      n3_loc = n3/np3;
      if (n1_loc==n2_loc&&n2_loc==n3_loc) break;
      if (n1/np1>=n2/np2) 
      {
        np1++;
      }
      else
      {
        np2++;
      }
      np3=np/(np1*np2);
    }
    if (np1*np2*np3!=np || 
        n1_loc!=n2_loc || n2_loc!=n3_loc ||
        np1*n1_loc!=n1 || np2*n2_loc!=n2 || np3*n3_loc!=n3 )
    {
        PHIST_SOUT(PHIST_WARNING,"our simple cartesian partitioning is not\n"
                                 "sophisticated enough for this setup, reverting\n"
                                 "to standard linear partitioing\n");
        n1=-2;
    }
    else
    {
      int rem=pid;
      pid1=rem%np1;
      rem=(rem-pid1)/np1;
      pid2=rem%np2;
      rem=(rem-pid2)/np2;
      pid3=rem%np3;
      global_offset = pid*(n1_loc*n2_loc);
      off1=pid1*n1_loc;
      off2=pid2*n2_loc;
      off3=pid3*n3_loc;
      PHIST_SOUT(PHIST_VERBOSE,"3D partitioning: %d x %d x %d domains of size %d x %d %d\n",
                np1, np2, np3, n1_loc, n2_loc, n3_loc);
    }
  }
  else if (n1==-2)
  {
    // partitioning failed
    return row;
  }
  else if (arg2==-1)
  {
    // iperm: convert new GID to original one
    ghost_gidx_t rem=row-global_offset;
    int i1=rem%n1_loc;
    rem=(rem-i1)/n1_loc;
    int i2=rem%n2_loc;
    rem=(rem-i2)/n2_loc;
    int i3=rem%n3_loc;
  
    return ((off3+i3)*n2_loc+(off2+i2))*n1_loc+i1;
  }
  else if (arg2==+1)
  {
    // perm: convert old GID to new one
    ghost_gidx_t rem=row;
    // global indices (i1,i2,i3)
    int i1=rem%n1;
    rem=(rem-i1)/n1;
    int i2=rem%n2;
    int i3=(rem-i2)/n2;
    // local indices (j1,j2,j3) on new partition
    int j1=i1%n1_loc;
    int j2=i2%n2_loc;
    int j3=i3%n3_loc;
    // on which partition is this node now?
    int p1=i1/n1_loc;
    int p2=i2/n2_loc;
    int p3=i3/n3_loc;
    int p = (p3*np2+p2)*np1+p1;
    ghost_gidx_t offset = p*(ghost_gidx_t)(n1_loc*n2_loc*n3_loc); 
    return offset + ((j3*n2_loc)+j2)*n1_loc+j1;
  }
  return -1; // shouldn't get here
}

