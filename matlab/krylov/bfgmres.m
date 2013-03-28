function [X,flag,relres,iter,resvec]=bfgmres( A,B,X0,opts,M)
%                                                        
% function X=bfgmres( A,B,opts,M)                       
%                                                        
% Right-preconditioned flexible block GMRES(m).          
% Implementation follows Saad (Alg. 9.6, p273),          
% with double CGS orthogonalization, adapted for multi-  
% ple RHS.                                               
%                                                        
% args                                                   
% ~~~~                                                   
%                                                        
% A and M may be anything accepted by apply_op.          
%                                                        
% options                                                
% ~~~~~~~                                                
%                                                        
% opts.tol - convergence tolerance                       
% opts.m - max size of Krylov subspace before restart.   
% opts.maxIter - maximum total number of iterations      
% opts.relax   - if true, orthogonality of the basis     
%                vectors is seen in a relaxed way, for   
%                instance the orthog() function will use 
%                only a single classical Gram-Schmidt    
%                step (default is false).                
%                                                        

tol=getopt(opts,'tol',1e-8);
maxIter=getopt(opts,'maxIter',300);
m=getopt(opts,'m',maxIter);
m=min(m,maxIter);
debug=getopt(opts,'debug',false);
relax=getopt(opts,'relax',false);

flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

n=size(B,1);
k=size(B,2);
V=zeros(n,m*k);
Z=V;

% stores the upper (H)essenberg matrix for each of the
% k systems solved, H is transformed 'on-the-fly' to  
% upper triangular matrix R by Givens-rotations.
H=zeros(k*(m+1),k*m);

% Givens or Householder coefficients for the
% QR-decomposition of H
trans=zeros(2,m,k);

bnorm = bvnorm(B,2);

% returns the column indices for block j in V
idx = @(i) [(i-1)*k+1:i*k];

% initial guess
if (isempty(X0))
  X = zeros(n,k);
else
  X=X0;
end

while (1)
  R=B-apply_op(X,A);
  rs=zeros((m+1)*k,k);
  beta = bvnorm(R,2);

  if (debug)
    disp(['(re-)start: actual residual(s): ',num2str(beta./bnorm)]);
  end
  if (max(beta/bnorm)<=tol)
    flag=0;
    break;
  end
  if (iter>maxIter)
    flag=-1;
    disp('max iter exceeded');
    break;
  end
  [V0,rs0]=orthog(R);
  V(:,idx(1))=V0;
  rs(1:k,1:k)=rs0;
  % Arnoldi - build orthogonal basis V and upper Hessenberg matrix H
  for j=1:m
    iter = iter + 1;
    jdx=idx(j);
    Z(:,idx(j))=apply_op(V(:,idx(j)),M);
    W=apply_op(Z(:,idx(j)),A);
    % two steps of CGS
    [Vcol,Hcol]=orthog(V(:,1:j*k),W,relax);
    V(:,idx(j+1)) = Vcol;
    [Hcol,transj,rsjjp1]=updateQR(H(1:j*k,1:(j-1)*k),...
        Hcol,trans(:,1:j-1,:),rs([idx(j),idx(j+1)],:));
    H(1:k*(j+1),jdx)=Hcol;
    trans(:,j,:)=transj;
    rs([idx(j),idx(j+1)],:)=rsjjp1;
    if (debug)
      disp('transformed H');
      H(1:k*(j+1),1:j*k);
      disp('rs');
      rs(1:(j+1)*k,:)
    end
    relres = abs(diag(rs(idx(j+1),:))).'./bnorm;
    resvec=[resvec;relres];
    print_iter(iter,relres);
    if (max(relres)<=tol)
      %disp('convergence - exit Arnoldi');
      break;
    end
  end % Arnoldi
  % H is upper triangular. We use triu here because
  % in the final implementation we might want to store
  % the Householder transformations in tril(H).
  y=triu(H(1:j*k,1:j*k))\rs(1:j*k,:);
  X = X+Z(:,1:j*k)*y;
end % while


end %function gmres
