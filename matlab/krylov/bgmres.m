function [X,flag,relres,iter,resvec]=bgmres( A,B,X0,opts,M)
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
% The block size is determined by size(B,2) right now.   
%

flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

n=size(B,1);
k=size(B,2);

idx = @(j) ((j-1)*k+1):(j*k);

tol=getopt(opts,'tol',1e-8);
maxIter=getopt(opts,'maxIter',300);
m=getopt(opts,'m',maxIter);
m=min(m,maxIter);
debug=getopt(opts,'debug',false);
relax=getopt(opts,'relax',false);

if (debug)
  tol
  maxIter
  m
end

V=zeros(n,m*k);

% stores the upper (H)essenberg matrix for each of the
% k systems solved, H is transformed 'on-the-fly' to  
% upper triangular matrix R by Givens-rotations.
H=zeros(k*(m+1),k*m);

% cos- and sin-terms for the Givens rotation
% to transform H to upper triangular form
%    1            c1 s1                  
%      c2 s2     -s1 c1                  
%..*  -s2 c2    *      1                 
%       I   1            1               
%                                        
% in the case k>1 we do not use          
% Givens but Householder transforms,     
% but they are stored in the same arrays 
%
% k=1: trans(1,:) = c, trans(2,:)=s.            
% k>1: trans(1,:) = beta, trans(2,:)=v(1),      
%      H(j,j)=alpha
%      H(j+[1:k],j) = v(2:end)
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
%    disp('max iter exceeded');
    break;
  end
  % orthogonalize initial residual
  [V0,rs0]=orthog(R);
  V(:,idx(1))=V0;
  rs(1:k,1:k)=rs0;
  % Arnoldi - build orthogonal basis V and upper Hessenberg matrix H
  for j=1:m
    iter = iter + 1;
    jdx=idx(j);
    if (debug)
      j
      jdx
    end
    W=apply_op(apply_op(V(:,jdx),M),A);
    % Arnoldi update
    [Vcol,Hcol]=orthog(V(:,1:j*k),W,relax);
    V(:,idx(j+1))=Vcol;
    % update QR factorization of H
    [Hcol,transj,rsjjp1]=updateQR(H(1:k*j,1:k*(j-1)),Hcol,...
        trans(:,1:j-1,:),rs([idx(j),idx(j+1)],:));
    H(1:k*(j+1),jdx)=Hcol;
    trans(:,j,1:k)=transj;
    rs([idx(j),idx(j+1)],:)=rsjjp1;
    
    if (debug)
        disp('transformed H');
        H(1:k*(j+1),1:j*k)
        disp('rs');
      rs(1:(j+1)*k,:)
    end
    relres = abs(diag(rs(idx(j+1),:))).'./bnorm;
    resvec=[resvec;relres];
    print_iter(iter,relres);    
    if (max(relres)<=tol || iter>=maxIter)
      %disp('convergence - exit Arnoldi');
      break;
    end
  end % Arnoldi
  y = triu(H(1:k*j,1:k*j))\rs(1:k*j,:);
  X = X+apply_op(V(:,1:j*k)*y,M);
end % while

end %function gmres
