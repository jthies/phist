function [X,flag,relres,iter,resvec]=bdqgmres( A,B,X0,opts,M)
%                                                        
% function X=bdqgmres( A,B,opts,M)                   
%                                                        
% Right-preconditioned block DQGMRES.                    
%                                                        
% The Direct Quasi-GMRES(m) method does not restart but  
% truncates the orthogonlization, i.e. we orthogonalize  
% only to the last m vectors in the basis. DQGMRES is an 
% alternative to F(lexible)GMRES as its inexact orthogo- 
% nalization allows a varying preconditioner. It is also 
% attractive because the Gram-Schmidt process has the    
% same length in each iteration j (for j>=m), which faci-
% litates vectorization/optimization.                    
%                                                        
% See Saad, Alg. 6.13 on p. 174.                         
%                                                        
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
% opts.m - number of vectors in basis (truncation param) 
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
debug=getopt(opts,'debug',false);
relax=getopt(opts,'relax',false);

flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

n=size(B,1);
k=size(B,2);

V=zeros(n,(maxIter+1)*k);
P=zeros(n,(maxIter+1)*k);

% stores the upper (H)essenberg matrix for each of the
% k systems solved, H is transformed 'on-the-fly' to  
% upper triangular matrix R by Givens-rotations.
H=zeros(k*(maxIter+2),k*(maxIter+1));

idx = @(jj) (jj-1)*k+1:jj*k; % as in bgmres, to index MV jj in a block like V

% indexing functions to get the location of V(:,j-m+1:j).
% We periodically overwrite the memory, so this is a bit tricky.

% for debugging - do not overwrite memory
%transIDX=@(jj) jj;
%ind_j = @(jj) max(jj-m+1,1):jj;
%ind_jk= @(jj,kk) kron((ind_j(jj)-1)*k,ones(1,length(kk)))+...
% kron(ones(1,length(ind_j(jj))),kk);

% periodically overwrite memory as we only 
% orthogonalize against the last m vectors
 transIDX=@(jj) mod(jj-1,m+2)+1;
 ind_j = @(jj) transIDX([max(jj-m+1,1):jj]);
 ind_jk= @(jj,kk) kron((ind_j(jj)-1)*k,ones(1,length(kk)))+...
                 kron(ones(1,length(ind_j(jj))),kk);
                                

% coefficients for Givens or Householder transforms of H
trans=zeros(2,maxIter,k);
rs=zeros((maxIter+1)*k,k);

bnorm = bvnorm(B,2);

% initial guess
if (isempty(X0))
  X = zeros(n,k);
else
  X=X0;
end

  R=B-apply_op(X,A);
  [V0,rs0]=orthog(R);
  rs(1:k,1:k)=rs0;
  V(:,idx(1))=V0;

  % Arnoldi - build (semi-)orthogonal basis V and upper Hessenberg matrix H
  for iter=1:maxIter
    j=transIDX(iter);
    jrange = ind_j(iter); % replaces 1:j or j-m+1:j for iter>m
    jkrange = ind_jk(iter,1:k); % replaces (1:k,j-m+1:j)
    jp1=transIDX(iter+1);
    jm1range = ind_j(iter-1); % replaces j-m:j-1
    jm1krange= ind_jk(iter-1,1:k);
    jp1range = [jrange,jp1]; % replaces j-m+1:j+1
    jp1krange= [jkrange,idx(jp1)];
    jkrangeX= [jm1krange,idx(j),idx(jp1)]; % j-m:j+1
    H(idx(jp1),:)=0;
    W=apply_op(apply_op(V(:,idx(j)),M),A);
    % orthogonalize
    [Vcol,p_Hcol]=orthog(V(:,jkrange),W,relax);
    V(:,idx(jp1))=Vcol;
    H(jp1krange,idx(j))=p_Hcol;
    rs(idx(jp1),:)=0;
    % update QR-decomposition of H
    [p_Hcol,transj,rsjjp1]=updateQR( H([jm1krange,idx(j)],jm1krange),...
                                    H(jkrangeX,idx(j)),...
                                    trans(:,jm1range,:),...
                                    rs([idx(j),idx(jp1)],:));
    trans(:,j,:)=transj;
    rs([idx(j),idx(jp1)],:)=rsjjp1;
    H(jkrangeX,idx(j))=p_Hcol;
    if (debug)
      disp('transformed H');
      H(jp1krange,jkrange)
      disp('rs');
      rs(jp1krange,:)
    end
    
    % add search vector P
    %Q=P(:,jm1krange)*H(jm1krange,idx(j));
    %Q=zeros(n,k);
    %for jj=jm1range
    %  Q=Q+P(:,idx(jj))*H(idx(jj),idx(j));
    %end
    %P(:,idx(j)) = V(:,idx(j)) - Q/triu(H(idx(j),idx(j)));
    P(:,idx(j)) = (V(:,idx(j)) - P(:,jm1krange)*H(jm1krange,idx(j)))/triu(H(idx(j),idx(j)));
    X=X+P(:,idx(j))*rs(idx(j),:);
    if (debug)
      disp('update');
      disp('V rows 1:4');
      V(1:4,jm1krange)
      disp('P rows 1:4')
      P(1:4,jm1krange)
      disp('X rows 1:4')
      X(1:4,:)
    end
    relres = abs(diag(rs(idx(jp1),:))).'./bnorm;
    resvec=[resvec;relres];
    print_iter(iter,relres);
%    disp(['    ',num2str(bvnorm(A*X-B)./bnorm)]);
    if (max(relres)<=tol)
      break;
    end
  end % Arnoldi

end %function gmres
