
void SUBR(pfgmres_init_args)(TYPE(pfgmres_args)* args,
        TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) rprec_op)
  {
  
  }      

// this variant of GMRES is intended for the parallel
// solution of multiple linear systems of the form
// (A-sigma_j I)X_j=B_j, j=1..num_sys. The number of
// vectors in X and B must be the same and a multiple
// of num_sys (B has num_sys*nb vectors). The rhs for 
// system j is assumed to be in B(:,j*bs:(j+1)*bs-1)
void SUBR(pfgmres)(TYPE(pfgmres_args)* args)
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
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(b,&k,ierr),*ierr);

  // only implemented for 1 rhs per shift right now:
  PHIST_CHK_IERR(*ierr=1-k,*ierr);

  mvec_ptr_t vj,zj; // views
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

// create ghost tasks for the operation y = alpha*x + beta*y
mvec_add_mvec_arg_t axpy_args;
axpy_args.ierr=ierr;
ghost_task_t* axpy_task = ghost_task_init(GHOST_TASK_FILL_ALL, GHOST_TASK_LD_ANY,
        (void*(*)(void*))SUBR(mvec_add_mvec_q), (void*)&axpy_args,GHOST_TASK_DEFAULT);

//TODO - sometimes we could put multiple operations in the buffer before
//       posting the wait, but the current taskBuf implementation does not
//       support that

// compute the norm of the RHS vector B
MT bnorm;
// put mvec_norm2 in the buffer
PHIST_CHK_IERR(taskBuf_add(b,(void*)&bnorm,args->t_id, args->op_NRM2,ierr),*ierr);
PHIST_CHK_IERR(taskBuf_wait(taskBuf,args->t_id,ierr),*ierr);
// initial guess TODO - give choice to user. This could also be put in the queue
// as it does not involve communication
PHIST_CHK_IERR(taskBuf_add(NULL,x,args->t_id, args->op_RNDX,ierr),*ierr);
PHIST_CHK_IERR(taskBuf_wait(taskBuf,args->t_id,ierr),*ierr);

while (1)
  {

  ////  r= b - (A+shift*I)x
  // r = A*x
  PHIST_CHK_IERR(taskBuf_add(args->taskBuf,x,r,args->t_id,args->op_AX,ierr),*ierr);
  PHIST_CHK_IERR(taskBuf_wait(taskBuf,args->t_id,ierr),*ierr);

  // r=r+shift*x via queue
  axpy_args.alpha=args->shift;
  axpy_args.beta=st::one();
  axpy_args.X=x;
  axpy_args.Y=r;
  ghost_task_add(axpy_task);
  ghost_task_wait(axpy_task);
  
  // put r=b-r in the queue
  axpy_args.alpha=st::one();
  axpy_args.beta=-st::one();
  axpy_args.X=b;
  axpy_args.Y=r;
  ghost_task_add(axpy_task);
  ghost_task_wait(axpy_task);
  
  for (int i=0; i<m+1;i++) rs[i]=st::zero();

  PHIST_CHK_IERR(taskBuf_add(r,(void*)&beta,args->t_id, args->op_NRM2,ierr),*ierr);
  PHIST_CHK_IERR(taskBuf_wait(taskBuf,args->t_id,ierr),*ierr);

  PHIST_OUT(2,"(re-)start: actual residual: %8.4g\n",beta/bnorm);

  if (beta/bnorm<=tol)
    {
    *ierr=0;
    break;
    }
  if (iter>args->max_iter)
    {
    *ierr=-1;
    PHIST_MSG(1,"max iter exceeded");
    break;
    }
  // get a view of the first column of V
  PHIST_CHK_IERR(SUBR(mvec_view_block)(V,vj,0,k,ierr),*ierr);

  // normalize r_0 and assign it to v_0
  //TODO - everything has to go in the queue...
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V,r,0,k,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_normalize)(vj,&rs[0],ierr),*ierr);
  
  // Arnoldi - build orthogonal basis V and upper Hessenberg matrix H
  for (int j=0; j<m; j++)
    {
    args->num_iter++;

    // get vector views
    PHIST_CHK_IERR(SUBR(mvec_view_block)(V,vj,0,k,ierr),*ierr);
    PHIST_CHK_IERR(SUBR(mvec_view_block)(Z,zj,0,k,ierr),*ierr);

    // TODO - better overlapping of communication and computation
    //{
    // z_j = M\v_j
    PHIST_CHK_IERR(taskBuf_add(vj,zj,args->t_id, args->op_RPRECX,ierr),*ierr);
PHIST_CHK_IERR(taskBuf_wait(taskBuf,args->t_id,ierr),*ierr);
    // w = op*z_j
    PHIST_CHK_IERR(taskBuf_add(zj,w,args->t_id, args->op_AX,ierr),*ierr);
    PHIST_CHK_IERR(taskBuf_wait(taskBuf,args->t_id,ierr),*ierr);
    
    // w = w + shift*z_j in queue
    axpy_args.alpha=args->shift;
    axpy_args.beta=st::one();
    axpy_args.X=zj;
    axpy_args.Y=w;
    ghost_task_add(axpy_task);
    ghost_task_wait(axpy_task);

  // TODO - maybe orthog should use the queue. That way
  // any kernel operations are done through the queue, 
  // but not algorithms like Gram-Schmidt.
  
    //[Vcol,Hcol]=orthog(V(:,1:j*k),W,relax);
    //V(:,idx(j+1)) = Vcol;
    //[Hcol,transj,rsjjp1]=updateQR(H(1:j*k,1:(j-1)*k),...
    //    Hcol,trans(:,1:j-1,:),rs([idx(j),idx(j+1)],:));
    //H(1:k*(j+1),jdx)=Hcol;
    //trans(:,j,:)=transj;
    //rs[j,j+1],:)=rsjjp1;
    relres = std::abs(rs[j+1])/bnorm;
    PHIST_OUT(0,"%d %8.4g\n",arg->num_iters,relres);
    if (relres<=tol)
      {
      PHIST_OUT(3,"convergence - exit Arnoldi");
      break;
      }
    } // Arnoldi
  // H is upper triangular. We use triu here because
  // in the final implementation we might want to store
  // the Householder transformations in tril(H).
  //TODO - lapack routine xtrsm
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
