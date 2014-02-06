function [x,flag,relres,iter,resvec]=carp_cg( A,b,x0,opts)
%                                                        
% function X=bgmres( A,B,X0,opts,M)                      
%                                                        
% Right-preconditioned block GMRES(m).                   
% Implementation follows Saad (Alg. 9.5, p270),          
% with single Block-CGS orthogonalization, adapted for   
% multiple RHS.                                          
%                                                        
% args                                                   
% ~~~~                                                   
%                                                        
% A should be a sparse matrix for now
%                                                        
% options                                                
% ~~~~~~~                                                
%                                                        
% opts.tol - convergence tolerance                       
% opts.maxIter - max number of iterations
% opts.omega - relaxation parameter
% opts.Precond - preconditioning operator
%

verbose=getopt(opts,'verbose',true);

itprint=10;
itcheck=1;

flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

n=size(A,1);

tol=getopt(opts,'tol',1e-8);
maxIter=getopt(opts,'maxIter',300);
omega=getopt(opts,'omega',1.7);
M=getopt(opts,'Precond',speye(n));

fprintf('CARP-CG tol: %4.2e\n',tol);

debug=getopt(opts,'debug',false);


if (debug)
  tol
  maxIter
end

nrm_b=norm(b);
reltol2=tol*tol*nrm_b*nrm_b;

nrms_ai2=zeros(n,1);
for i=1:n
  nrms_ai2(i) = A(i,:)*A(i,:)';
end

x=x0;
r=dkswp(A,b,x,omega,nrms_ai2)-x;
z=apply_op(r,M);
p=z;

r2_new = r'*z;

bnul=zeros(n,1);

disp(sprintf('%d\t%e\t%e',0,sqrt(r2_new),sqrt(r2_new)/nrm_b));
for k=1:maxIter
  q=p-dkswp(A,bnul,p,omega,nrms_ai2);
  alpha = (r'*z)/(p'*q);
  x=x+alpha*p;
  if (mod(k-1,itcheck)==0)
    nrm_r = norm(A*x-b);
    if (mod(k-1,itprint)==0)
      disp(sprintf('%d\t%e\t%e',k,nrm_r,nrm_r/nrm_b));
    end
    relres=nrm_r/nrm_b;
  %  resvec=[resvec,nrm_r];
    if (nrm_r<tol)
      break;
    end
  end
  r=r-alpha*q;
  z=apply_op(r,M);
  r2_old=r2_new;
  r2_new=r'*z;
  beta=r2_new/r2_old;
  p=z+beta*p;
  relres=sqrt(r2_old)/nrm_b;
  %fprintf('\t%d\t%e\n',k,relres);
  %if (r2_old<reltol2) 
  %  break;
  %end
end
iter=k;

if iter>=maxIter
  flag=1;
else
  flag=0;
end

end
