#include "matfuncs.h"

ghost_gidx_t iperm2d( ghost_gidx_t row, gidx_t arg2 )
{
  static int n1=-1, n2;
  static int np, np1, np2, pid, pid1, pid2;
  static int n1_loc, n2_loc;
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
      pid1 = np%np1;
      pid2 = (np-pid1)/np1;
      global_offset = pid*(n1_loc*n2_loc);
      PHIST_SOUT(PHIST_VERBOSE,"local partition size: %dx%d\n",
        n1_loc, n2_loc);
    }
  
  }
  else if (n1==-2)
  {
    // partitioning failed
    return row;
  }
  else
  {
    int i1 = row%n1;
    int i2 = (row-i1)/n1;
    int j1=i1%n1_loc; // local indices
    int j2=i2%n2_loc;
    return global_offset + j2*n1_loc+j1;
  }
  
}


ghost_gidx_t iperm3d( ghost_gidx_t row, gidx_t arg2, gidx_t arg3 )
{
  static int n1=-1, n2, n3;
  static int np, np1, np2, np3, pid, pid1, pid2, pid3;
  static int n1_loc, n2_loc, n3_loc;
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
    while (np1<=np && n2<=np)
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
      global_offset = pid*(n1_loc*n2_loc);
      PHIST_SOUT(PHIST_VERBOSE,"local partition size: %dx%dx%d\n",
        n1_loc, n2_loc, n3_loc);
    }
  
  }
  else if (n1==-2)
  {
    // partitioning failed
    return row;
  }
  else
  {
    ghost_gidx_t rem=row;
    int i1=rem%n1;
    rem=(rem-i1)/n1;
    int i2=rem%n2;
    rem=(rem-i2)/n2;
    int i3=rem%n3;
  
    int j1=i1%n1_loc; // local indices
    int j2=i2%n2_loc;
    int j3=i3%n3_loc;
    return global_offset + (j3*n2_loc+j2)*n1_loc+j1;
  }
  
}

