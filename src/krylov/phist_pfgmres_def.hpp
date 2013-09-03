// this variant of GMRES is intended for the parallel
// solution of multiple linear systems of the form
// (A-sigma_j I)X_j=B_j, j=1..num_sys. The number of
// vectors in X and B must be the same and a multiple
// of num_sys (B has num_sys*nb vectors). The rhs for 
// system j is assumed to be in B(:,j*bs:(j+1)*bs-1)
void SUBR(pbgmres)(TYPE(pfgmres_args)* args)
  {
#include "phist_std_typedefs.hpp"  
  int *ierr = &args->ierr;
  *ierr=0;

  int m = args->num_blocks;

  MT relres=1.0;
  args->num_iter = 0;
  mvec_ptr_t x = args->lhs;
  const_mvec_ptr_t b = args->rhs;

  int k;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(B,&k,ierr),*ierr);

  // only implemented for 1 rhs per shift right now:
  PHIST_CHK_IERR(*ierr=1-k,*ierr);

  mvec_ptr_t v,z; // views
  mvec_ptr_t V, Z;
  mvec_ptr_t r, w;
  const_map_ptr_t map;
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(x,&map,ierr),*ierr);
  PHIST_CHK_IERR(map_get_comm(map,&comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&r,map,k,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&w,map,k,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&V,map,m*k,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_create)(&Z,map,m*k,ierr),*ierr);
  // stores the upper (H)essenberg matrix for each of the
  // k systems solved, H is transformed 'on-the-fly' to  
  // upper triangular matrix R by Givens-rotations.
  sdMat_ptr_t H;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&H,m,m,comm,ierr),*ierr);

// cos- and sin-terms for the Givens rotation
// to transform H to upper triangular form
//    1            c1 s1                  
//      c2 s2     -s1 c1                  
//..*  -s2 c2    *      1                 
//       I   1            1               
//                                        
//
ST givS=new ST[m];
ST givC=new ST[m];
ST rs = new ST[m+1];

MT beta;

// compute the norm of the RHS vector B
MT bnorm;
// put mvec_norm2 into queue
{

}
// initial guess TODO - give choice to user
// TODO - put random() into the queue

while (1)
  {

  //  R=B-op*X
  PHIST_CHK_IERR(taskBuf_add(args->taskBuf,X,R,args->t_id,args->op_OPX,ierr),*ierr);
  // TODO put R=B-R in the queue
  for (int i=0; i<m+1;i++) rs[i]=st::zero();
  // TODO - put this in the queue
  beta = bvnorm(R,2);

  if (debug)
    {
    disp(['(re-)start: actual residual(s): ',num2str(beta./bnorm)]);
    }
  if (max(beta/bnorm)<=tol)
    {
    *ierr=0;
    break;
    }
  if (iter>maxIter)
    {
    *ierr=-1;
    PHIST_MSG(1,"max iter exceeded");
    break;
    }
  // TODO - orthog should use the queue
  //[V0,rs0]=orthog(R);
  //V(:,0)=V0;
  //rs[0]=rs0;
  // Arnoldi - build orthogonal basis V and upper Hessenberg matrix H
  for (int j=0; j<m; j++)
    {
    args->num_iter++;
    // TODO - better overlapping of communication and computation
    Z(:,idx(j))=apply_op(V(:,idx(j)),M);
    W=apply_op(Z(:,idx(j)),A);
    // two steps of CGS
    [Vcol,Hcol]=orthog(V(:,1:j*k),W,relax);
    V(:,idx(j+1)) = Vcol;
    [Hcol,transj,rsjjp1]=updateQR(H(1:j*k,1:(j-1)*k),...
        Hcol,trans(:,1:j-1,:),rs([idx(j),idx(j+1)],:));
    H(1:k*(j+1),jdx)=Hcol;
    trans(:,j,:)=transj;
    rs([idx(j),idx(j+1)],:)=rsjjp1;
    relres = abs(diag(rs(idx(j+1),:))).'./bnorm;
    PHIST_OUT(0,"%d %8.4g\n",arg->num_iters,relres);
    if (max(relres)<=tol)
      {
      PHIST_OUT(3,"convergence - exit Arnoldi");
      break;
      }
    } // Arnoldi
  // H is upper triangular. We use triu here because
  // in the final implementation we might want to store
  // the Householder transformations in tril(H).
  //y=triu(H(1:j*k,1:j*k))\rs(1:j*k,:);
  //X = X+Z(:,1:j*k)*y;


  } //end while

delete [] rs;
delete [] givS;
delete [] givC;

PHIST_CHK_IERR(SUBR(mvec_delete)(r,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_delete)(w,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_delete)(V,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_delete)(Z,ierr),*ierr);
PHIST_CHK_IERR(SUBR(sdMat_delete)(H,ierr),*ierr);

}

function [Hcol,transj,rsjjp1]=updateQR(Hj,Hcol,trans,rsjjp1)
//                                                               
//function [Hcol,transj,rsjjp1]=updateQR(Hj,Hcol,trans,rsjjp1)   
//                                                               
// input: previous j x j-1 (j*k x (j-1)*k) tansformed            
//       H-matrix Hj,                                            
//        new column Hcol (dim. (j+1)*k x k)                     
//        transformation coefficients (Givens or HH) for         
//        previous j-1 rows of H                                 
//        right-hand side to be transformed along with H [j,j+1] 
//                                                        
// output Hcol and rs overwrite inputs, trans(:,j)=transj.
//                                                        

j=size(Hj,1);

// apply previous (j-1) transformations to columns j
  for jj=1:size(trans,2) % apply Givens rotation
    {
    htmp = trans(1,jj)*Hcol(jj) + ...
           trans(2,jj)*Hcol(jj+1);
    Hcol(jj+1) = -conj(trans(2,jj))*Hcol(jj) + trans(1,jj)*Hcol(jj+1);
    Hcol(jj)   = htmp;
    }

// update QR factorization of H

  // new Givens rotation for eliminating H(j+1,j)
  [cs,sn,htmp]=give(Hcol(j),Hcol(j+1));
  transj(1,1)=cs;
  transj(2,1)=sn;
  // eliminate H(j+1,j)
  Hcol(j) = htmp;
  Hcol(j+1) = 0.0;
  // apply to RHS
  tmp=cs*rsjjp1(idx(1));
  rsjjp1(idx(2)) = -conj(sn).*rsjjp1(idx(1));
  rsjjp1(idx(1))=tmp;

} // updateQR
