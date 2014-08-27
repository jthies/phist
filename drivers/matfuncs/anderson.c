#include "matfuncs.h"

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

int anderson( ghost_idx_t row, ghost_idx_t *nnz, ghost_idx_t *cols, void *vals){

	static ghost_idx_t nx = 4 ;
	ghost_idx_t ny=nx;
	ghost_idx_t nz=nx;
	ghost_idx_t N = nx*ny*nz;

	ghost_idx_t           max_row_nnz  = 7;


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

		return 0;
	}else if ( row == -2) 
	{
		nx = nnz[0];
		ny = nx;
		nz = nx;

		//printf("W,L = %d, %d\n",W,L);
		//printf("nnz = %d\n",max_row_nnz);
		return 0;
	}

	printf("################### error \n");
	*nnz = -1; //  error
	return 1;
}


