function [V,H]=arnoldi(A,v0,opts)
% function [V,H]=arnoldi(A,v0,opts)                             
% (Block-)Arnoldi decomposition, A*V(:,1:m*k)=V(:,1:(m+1)*k)*H  
% with V'*V=I, H (block-)upper Hessenberg. The block size is de-
% termined by the number of columns in the starting vector v0.  
% v0 need not be orthonormal when passed in.                    
%                                                               
% Up to now this function does not compute any eigenvalues or   
% such. Valid options are:                                      
%                                                               
% opts.maxIter - number of steps to take
%                                                               


n=size(v0,1);
k=size(v0,2);

maxIter=getopt(opts,'maxIter',30);
m=min(floor(n/k)-1,maxIter);

idx=@(j) ((j-1)*k+1):(j*k);

% generate k orthogonal starting vectors
V = orthog(v0);

H=zeros((m+1)*k,m*k);

for j=1:m
  W=apply_op(V(:,idx(j)),A);

  % orthogonalize
  [Vnew,hcol]=orthog(V,W);
  V=[V, Vnew];
  H(1:max(idx(j+1)),idx(j))=hcol;
end

end
