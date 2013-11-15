function [X,flag,relres,iter,resvec]=simple_gmres( A,B,X0,opts,M)
% single vector variant of bgmres.m for debugging the phist implementation
% in src/jada/phist_simple_gmres*

verbose=getopt(opts,'verbose',true);

flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

n=size(B,1);

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

V=zeros(n,m);

% stores the upper (H)essenberg matrix for each of the
% k systems solved, H is transformed 'on-the-fly' to  
% upper triangular matrix R by Givens-rotations.
H=zeros(m+1,m);

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
trans=zeros(2,m);

bnorm = bvnorm(B,2);

% initial guess
if (isempty(X0))
  X = zeros(n,1);
else
  X=X0;
end

while (1)
  R=B-apply_op(X,A);
  rs=zeros(m+1,1);
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
  rs(1)=rs0;
  V(:,1)=V0;
  % Arnoldi - build orthogonal basis V and upper Hessenberg matrix H
  for j=1:m
    iter = iter + 1;
    if (debug)
      j
      j
    end
    W=apply_op(apply_op(V(:,j),M),A);
    % Arnoldi update
    [Vcol,Hcol]=orthog(V(:,1:j),W,relax);
    V(:,j+1)=Vcol;
    % update QR factorization of H

    % apply previous (j-1) transformations to columns j
    for jj=1:j-1 % apply Givens rotation
      htmp = trans(1,jj)*Hcol(jj) + ...
           trans(2,jj)*Hcol(jj+1);
      Hcol(jj+1) = -conj(trans(2,jj))*Hcol(jj) + trans(1,jj)*Hcol(jj+1);
      Hcol(jj)   = htmp;
    end

    % new Givens rotation for eliminating H(j+1,j)
    [cs,sn,htmp]=give(Hcol(j),Hcol(j+1));
    trans(1,j)=cs;
    trans(2,j)=sn;
    % eliminate H(j+1,j)
    Hcol(j) = htmp;
    Hcol(j+1)=0;
    % apply to RHS
    tmp=cs*rs(j);
    rs(j+1) = -conj(sn).*rs(j);
    rs(j)=tmp;

    H(1:j+1,j)=Hcol;
    
    if (debug)
        disp('transformed H');
        H(1:k*(j+1),1:j)
        disp('rs');
      rs(1:j+1,:)
    end
    relres = abs(rs(j+1))./bnorm;
    resvec=[resvec;relres];
    if (verbose)
      print_iter(iter,relres);
    end
    if (max(relres)<=tol || iter>=maxIter)
      %disp('convergence - exit Arnoldi');
      break;
    end
    %disp(sprintf('j=%d, cs=%8.4f sn=%8.4f rs=%8.4f',j,trans(1,j),trans(2,j),rs(j)));
  end % Arnoldi
  y = triu(H(1:j,1:j))\rs(1:j);
  X = X+apply_op(V(:,1:j)*y,M);
end % while

end %function gmres
