function [Q,R1,R2]=orthog(varargin)
%                                                               
% [Q,R]=orthog(W)                                               
%                                                               
% Makes the columns of W mutually orthogonal, W=Q*R, Q'Q=I and  
% throws out the null space of W if W does not have full rank.  
%                                                               
% [Q,R1,R2]=orthog(V,W [, relax])                               
%                                                               
% First  orthogonalizes the columns of W against all vectors    
% in V using classical Gram-Schmidt, then orthogonalizes W in-  
% ternally. If W does not have full rank, we throw out the      
% null space, so that Q spans the same space as W. The process  
% is repeated for numerical stability, the relax flag is        
% currently ignored. V is assumed to have orthonormal columns   
% already.                                                      
%                                                               
% returns Q*R1=W-V*R2 such that V'*Q=0, Q'Q=I.                  
%                                                               

relax=false;
if (nargin==1)
  W=varargin{1};
  V=zeros(size(W,1),0);
elseif (nargin>1)
  V=varargin{1};
  W=varargin{2};
elseif (nargin<1 || nargin>3)
  error('invalid number of arguments');
end
if (nargin==3)
  relax=varargin{3};
end

if isempty(V)
  V=zeros(size(W,1),0);
end


n=size(V,1);
m=size(V,2);
k=size(W,2);

R=zeros(m+k,k);

% orthogonalize against V
R1=V'*W;
Z=W-V*R1;

% orthogonalize Z and compute numerical rank
[Q,R2]=qr(Z);

rank = length(find(abs(diag(R2))>eps));
disp(['rank V after first GS pass: ',num2str(rank)]);

if (relax==false)
  %Z=bvscal(Z,1./bvnorm(Z,2));
  R2=V'*Z;
  Z=Z-V*R2;
  R1=R1+R2;
end

% orthogonalize Z
[Q,R2]=qr(Z,0);

R2=bvscal(R2,znrm);

R(1:m,1:k)=R1;
R(m+1:m+k,1:k)=R2;

end
