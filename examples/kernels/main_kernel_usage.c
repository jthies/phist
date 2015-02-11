#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_macros.h"
#include <stdio.h>

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h" 
#endif

// in order to use the driver utils header, we first have to
// specify which data type we will use (real double (D) here)
#include "phist_gen_d.h"
#include "phist_driver_utils.h"

/* This is a first example program to introduce 
   the PHIST concept of kernel operations independent
   of actual data structures.
   
   The program first reads in a sparse matrix from a file
   or generates it from a row function. Then it constructs
   two blocks of vectors with 4 columns each, fills one of them
   with ones and multiplies it by the matrix into
   the other.
   
   The process is then repeated with a "view" of the first two columns
   of the vectors, and then with the second two. In the second matrix-vector
   product, we subtract the target vectors (y<-A*x-y), so the end result should
   be 0, which we check by taking the col-wise 2-norm.
   */
int main(int argc, char** argv)
{
  // MPI info
  int rank, num_proc;
  // every PHIST function has as large argument "iflag", an integer return code that is 0 if 
  // the function was successful, negative if it failed, and positive if it issued a warning
  int iflag;

  // opague objects used in the underlying kernel lib
  
  // a wrapper for the MPI comm (may actually be something completely different if the 
  // kernel lib does not use MPI for communication)
  const_comm_ptr_t comm = NULL;
  // sparse matrix
  DsparseMat_ptr_t A;
  // defines the distribution of matrix/vector rows and columns
  const_map_ptr_t row_map,range_map,domain_map;
  // block vectors
  Dmvec_ptr_t x,y;
  // views of certain columns of x and y. These objects
  // can be used exactly like standard mvecs.
  Dmvec_ptr_t x_view, y_view;

  // we can "view" the raw data of vectors and manipulate it, but have to
  // be *very* careful to up- and download it to/from an accelerator if necessary
  double *x_val, *y_val;
  lidx_t nloc_x, nloc_y;
  int nvec_x,nvec_y;
  lidx_t lda_x, lda_y;

  int i,j;
  
  comm_ptr_t comm_world = NULL;
  double norms[4];

  // the PHIST_ICHK_IERR (PHIST_CHK_IERR in void functions) macros provide
  // some convenient checking of the return flag.
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  // this generates the default comm object, typically MPI_COMM_WORLD for MPI-based kernel libs
  PHIST_ICHK_IERR(phist_comm_create(&comm_world,&iflag),iflag);

  // print some usage info if the program is called without arguments
  if (argc<2)
  {
    // the PHIST_SOUT macro wraps printf on the root process and if the given level is at 
    // most PHIST_OUTLEV
    PHIST_SOUT(PHIST_ERROR,"Usage: %s <matrix>\n"
                           "       \"matrix\" is the matrix file to be read/created (compulsory)\n\n",
                                argv[0]);
        // we use this function below to read/generate the matrix, and calling it
        // with the "usage" string makes it print a useful usage message.
        phist_Dcreate_matrix(NULL,NULL,"usage",&iflag);
    return 1;
  }


  // this is a very convenient function for reading a matrix from a file or generating one 
  // of our scalable test cases (from ESSEX-PHYSICS or the builtin surrogate)
  PHIST_ICHK_IERR(phist_Dcreate_matrix(&A,comm_world,argv[1],&iflag),iflag);
    
  // the range map defines the data distribution and ordering of the LHS vectors 
  // (Y in a product Y=A*X)
  PHIST_ICHK_IERR(phist_DsparseMat_get_range_map(A, &range_map, &iflag),iflag);
  // the domain map defines the data distribution and ordering of the RHS vectors 
  // (X in a product Y=A*X)
  PHIST_ICHK_IERR(phist_DsparseMat_get_domain_map(A, &domain_map, &iflag),iflag);

  // this is the communicator actually used by the matrix object:
  PHIST_ICHK_IERR(phist_map_get_comm(range_map, &comm, &iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank, &iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc, &iflag),iflag);

  // now we know the data distribution, we can create our block vectors.
  // We store four vectors in each block x, y:
  PHIST_ICHK_IERR(phist_Dmvec_create(&x,domain_map,4,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_create(&y,range_map,4,&iflag),iflag);

  // find out the properties of the vectors:
  PHIST_ICHK_IERR(phist_Dmvec_my_length(x,&nloc_x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_my_length(y,&nloc_y,&iflag),iflag);

  PHIST_ICHK_IERR(phist_Dmvec_num_vectors(x,&nvec_x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_num_vectors(y,&nvec_y,&iflag),iflag);
  
  fprintf(stdout,"rank %d: x has local length %d and %d vectors\n",rank,(int)nloc_x,(int)nvec_x);
  //fprintf(stdout,"rank %d: y has local length %d and %d vectors\n",rank,(int)nloc_y,(int)nvec_y);
  
  // get pointers to the raw data. Beware that the data may be stored either in row- or 
  // column-major ordering in PHIST, depending on the macro PHIST_MVECS_ROW_MAJOR.
  PHIST_ICHK_IERR(phist_Dmvec_extract_view(x,&x_val,&lda_x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_extract_view(y,&y_val,&lda_y,&iflag),iflag);
  
  // set x=1
  PHIST_ICHK_IERR(phist_Dmvec_put_value(x,1.0,&iflag),iflag);

  // this is the "manual" way to set vector entries. For two
  // reasons it is very dangerous: if the host node is an accelerator,
  // you manually have to upload the new data after setting it. And if
  // the host node is e.g. an OpenMP node, you may be messing with the
  // NUMA placement of the vector entries.
#ifdef PHIST_MVECS_ROW_MAJOR
  for (i=0;i<nloc_y;i++)
    for (j=0;j<nvec_y;j++)
      y_val[i*lda_x+j]=1.0;
#else
  for (j=0;j<nvec;j++)
    for (i=0;i<nloc;i++)
      y_val[j*lda_x+i]=1.0;
#endif
  // after manually manipulating vector entries, we potentially need to upload them to
  // an accelerator. If not, this call will do nothing.
  PHIST_ICHK_IERR(phist_Dmvec_to_device(y,&iflag),iflag);

  // compute y=A*x
  PHIST_ICHK_IERR(phist_DsparseMat_times_mvec(1.0,A,x,0.0,y,&iflag),iflag);

  // compute and print 2-norms
  PHIST_ICHK_IERR(phist_Dmvec_norm2(y, norms,&iflag),iflag);

  PHIST_SOUT(PHIST_INFO,"after first MVM: cols contain A*ones.\n");
  for (i=0;i<4;i++)
  {
    PHIST_SOUT(PHIST_INFO,"%8.4g  ",norms[i]);
  }
  PHIST_SOUT(PHIST_INFO,"\n");
  
  // create views of the first two columns of x and y. It is *very* important
  // to set the result vectors to NULL if they are not yet allocated by some 
  // previous calls because PHIST otherwise assumes that you're trying to update
  // an existing view. The bounds given are inclusive, so we view the first two
  // vectors here.
  x_view=NULL;
  PHIST_ICHK_IERR(phist_Dmvec_view_block(x,&x_view,0,1,&iflag),iflag);
  y_view=NULL;
  PHIST_ICHK_IERR(phist_Dmvec_view_block(y,&y_view,0,1,&iflag),iflag);
  
  // now do the second matrix vector product on the first two columns,
  // subtracting the previous result
  PHIST_ICHK_IERR(phist_DsparseMat_times_mvec(1.0,A,x_view,-1.0,y_view,&iflag),iflag);

  // compute and print 2-norms
  PHIST_ICHK_IERR(phist_Dmvec_norm2(y, norms,&iflag),iflag);

  PHIST_SOUT(PHIST_INFO,"after second MVM: first cols now should be 0.\n");
  for (i=0;i<4;i++)
  {
    PHIST_SOUT(PHIST_INFO,"%8.4g  ",norms[i]);
  }
  PHIST_SOUT(PHIST_INFO,"\n");

  // update the views to view the second two columns. The kernel lib will rebuild or update
  // the mvecs, and no memory leak occurs.
  PHIST_ICHK_IERR(phist_Dmvec_view_block(x,&x_view,2,3,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_view_block(y,&y_view,2,3,&iflag),iflag);

  // again, subract A*x(:,2:3) from y(:,2:3)
  PHIST_ICHK_IERR(phist_DsparseMat_times_mvec(1.0,A,x_view,-1.0,y_view,&iflag),iflag);
  
  // compute and print 2-norms
  PHIST_ICHK_IERR(phist_Dmvec_norm2(y, norms,&iflag),iflag);
  
  PHIST_SOUT(PHIST_INFO,"after third MVM: all norms now should be 0.\n");
  for (i=0;i<4;i++)
  {
    PHIST_SOUT(PHIST_INFO,"%8.4g  ",norms[i]);
  }
  PHIST_SOUT(PHIST_INFO,"\n");
  
  // delete everything
  PHIST_ICHK_IERR(phist_Dmvec_delete(x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_delete(y,&iflag),iflag);
  PHIST_ICHK_IERR(phist_DsparseMat_delete(A,&iflag),iflag);

  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  }
