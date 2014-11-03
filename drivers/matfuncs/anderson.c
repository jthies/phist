#include "matfuncs.h"

#include <math.h>

#ifdef WRITE_MATRIX
#include <mpi.h>
#endif

#ifndef anderson_L
#define anderson_L 16.5
#endif

// generate a simple matrix representing a 3D 7-point stencil
// with random numbers on the diagonal (between -L/2 and L/2,
// where L=16.5 is fixed), -1 on the off-diagonals and periodic BC.

#ifdef MOD
#undef MOD
#endif

#define MOD(x,y) (((double)(y)==0.0)? (double)(x): ((double)(x) - floor((double)(x)/((double)(y)))*((double)(y))))


int gid2ijk(ghost_gidx_t gid, 
        int nx, int ny, int nz,
        int* i, int* j, int* k)
{
    if (gid<0 || gid>nx*ny*nz)
    {
      return -1;
    }
    ghost_gidx_t rem=gid;
    *i=MOD(rem,nx);
    rem=(rem-*i)/nx;
    *j=MOD(rem,ny);
    rem=(rem-*j)/ny;
    *k=MOD(rem,nz);
    return 0;
}

int ijk2gid(int nx, int ny, int nz,
        int i, int j, int k)
{
  int ii=i, jj=j, kk=k;
  
if (ii<0) ii+=nx;
  else if (ii>=nx) ii-=nx;

if (jj<0) jj+=ny;
  else if (jj>=ny) jj-=ny;

if (kk<0) kk+=nz;
  else if (kk>=nz) kk-=nz;
  
  return (kk*ny+jj)*nx+ii;
}

int anderson( ghost_gidx_t row, ghost_lidx_t *nnz, ghost_gidx_t *cols, void *vals){

	static ghost_lidx_t nx = 0 ;
	ghost_lidx_t ny=nx;
	ghost_lidx_t nz=nx;
	ghost_gidx_t N = nx*ny*nz;

	ghost_lidx_t           max_row_nnz  = 7;

#ifdef WRITE_MATRIX
        int me;
        MPI_Comm_rank(MPI_COMM_WORLD,&me);
        char fname[8];
        sprintf(fname,"mat%d.txt",me);
        FILE* deb=fopen(fname,"a");
#endif
        if ((row >-1 ) && (row <N)){       //  defined output -- write entries of #row in *cols and *vals
	                                   //                    return number of entries
                double * dvals = vals;

                int ii,jj,kk;
                gid2ijk(row,nx,ny,nz,&ii,&jj,&kk);

                double gamma=anderson_L;
                double V = gamma*( ((double)(rand()))/((double)(RAND_MAX)) -0.5);

                int i=0;
                dvals[i]=V; cols[i++]=row;
                dvals[i]=-1; cols[i++]=ijk2gid(nx,ny,nz,ii+1,jj,kk);
                dvals[i]=-1; cols[i++]=ijk2gid(nx,ny,nz,ii-1,jj,kk);
                dvals[i]=-1; cols[i++]=ijk2gid(nx,ny,nz,ii,jj+1,kk);
                dvals[i]=-1; cols[i++]=ijk2gid(nx,ny,nz,ii,jj-1,kk);
                dvals[i]=-1; cols[i++]=ijk2gid(nx,ny,nz,ii,jj,kk+1);
                dvals[i]=-1; cols[i++]=ijk2gid(nx,ny,nz,ii,jj,kk-1);
                *nnz = i;
#ifdef WRITE_MATRIX
        for (i=0;i<*nnz;i++)
        {
          fprintf(deb,"%d %d %24.16e\n",(int)row,(int)cols[i],dvals[i]);
          if (cols[i]<0 || cols[i]>=N)
          {
            fprintf(stdout,"proc %d, global row %d: column index %d found\n",me,row,cols[i]);
            fflush(stdout);
          }
        }
                
	fclose(deb);
#endif
        return 0;
        }

	// =================================================================
	else if ( row == -1) {
		matfuncs_info_t *info = vals;

		info->version   = 1;
		info->base      = 0;
		info->symmetry  = GHOST_SPARSEMAT_SYMM_GENERAL;
		info->datatype  = GHOST_DT_DOUBLE|GHOST_DT_REAL;
		info->nrows     = N;
		info->ncols     = N;
		info->row_nnz   = max_row_nnz;

		info->hermit    =  1;
		info->eig_info  =  0;
#ifdef WRITE_MATRIX
                fprintf(deb,"%ld %ld %ld\n",N,N,N*7);
   		fclose(deb);
#endif
		return 0;
	}else if ( row == -2) 
	{
		nx = nnz[0];
		ny = nx;
		nz = nx;
#ifdef WRITE_MATRIX
                fprintf(deb,"MatrixMarket matrix coordinate real general");
                fprintf(deb,"%Anderson model, %dx%dx%d\n",nx,ny,nz);
		fclose(deb);
#endif
		return 0;
	}

	printf("################### error \n");
	*nnz = -1; //  error
	return 1;
}


