function varargout=orthog(varargin)
%                                                               
% [Q,R]=orthog(W)                                               
%                                                               
% Makes the columns of W mutually orthogonal, W=Q*R, Q'Q=I and  
% throws out the null space of W if W does not have full rank.  
%                                                               
% [Q,Rcol]=orthog(V,W [, relax])                               
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
% returns Q*R1=W-V*R2 such that V'*Q=0, Q'Q=I. If only two out- 
% put args are given, Rcol=[R1;R2] is returned (useful for cre- 
% ating a Hessenberg-matrix in Arnoldi and such)                
%                                                               

relax=false;
if (nargin==1)
  W=varargin{1};
  V=zeros(size(W,1),0);
  if (nargout>2)
    error('usage: [Q,R]=orthog(W) or [Q,R1,R2]=orthog(V,W)');
  end
elseif (nargin>1)
  V=varargin{1};
  W=varargin{2};
  if (nargout>3)
    error('usage: [Q,R]=orthog(W) or [Q,R1,R2]=orthog(V,W)');
  end
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

% orthogonalize against V
R1=V'*W;
Z=W-V*R1;

% orthogonalize Z and compute numerical rank
[Q,R2]=qr(Z);

rank = length(find(abs(diag(R2))>eps));
%disp(['rank V after first GS pass: ',num2str(rank)]);

if (relax==false)
  R2=V'*Z;
  Z=Z-V*R2;
  R1=R1+R2;
end

% orthogonalize Z
[Q,R2]=qr(Z,0);

if (nargout>=1)
  varargout{1}=Q;
end
if (nargout==2)
  if (nargin==1)
    varargout{2}=R2;
  else
    varargout{2}=[R1;R2];
  end
else if (nargout==3)
  varargout{2}=R1;
  varargout{3}=R2;
end

end
