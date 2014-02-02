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

verbose=getopt(opts,'verbose',true);

itprint=10;

flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

n=size(A,1);

tol=getopt(opts,'tol',1e-8);
maxIter=getopt(opts,'maxIter',300);
omega=getopt(opts,'omega',1.7);
debug=getopt(opts,'debug',false);

if (debug)
  tol
  maxIter
end

nrm_b=norm(b);

nrms_ai2=zeros(n,1);
for i=1:n
  nrms_ai2(i) = A(i,:)*A(i,:)';
end

y=x0;
r=dkswp(A,b,y,omega,nrms_ai2)-y;
p=r;

r2_new = r'*r;

nul=zeros(n,1);

for k=1:maxIter
  q=p-dkswp(A,nul,p,omega,nrms_ai2);
  alpha = (r'*r)/(p'*q);
  y=y+alpha*p;
  if (mod(k,itprint)==0)
    nrm_r = norm(A*y-b);
    disp(sprintf('%d\t%e\t%e',k,nrm_r,nrm_r/nrm_b));
    relres=nrm_r/nrm_b;
    resvec=[resvec,nrm_r];
    if (nrm_r<tol*nrm_b)
      break;
    end
  end
  r=r-alpha*q;
  r2_old=r2_new;
  r2_new=r'*r;
  beta=r2_new/r2_old;
  p=r+beta*p;
  %disp(sprintf('%d\t%f',k,beta));
  if (beta<tol) 
    break;
  end
end
iter=k;
%x=A'*y;
x=y;
end
