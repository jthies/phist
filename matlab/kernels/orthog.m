function [Q,R]=orthog(varargin)
%                                                               
% [Q,R]=orthog(W)                                               
%                                                               
% Makes the columns of W mutually orthogonal, W=Q*R, Q'Q=I      
%                                                               
% [Q,R]=orthog(V,W [, relax])                                   
%                                                               
% First  orthogonalizes the columns of W against all vectors    
% in V using classical Gram-Schmidt. The GS process is repeated,
% unless relax=true is specified. V is assumed to have ortho-   
% normal columns already.                                       
% Then makes the columns of W mutually orthogonal (by a QR-     
% factorization),                                               
%                                                               
% returns Q*R=W such that V'*Q=0, Q'Q=I.                        
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

n=size(V,1);
m=size(V,2);
k=size(W,2);

R=zeros(m+k,k);

% orthogonalize against V
R1=V'*W;
Z=W-V*R1;
if (relax==false)
  %Z=bvscal(Z,1./bvnorm(Z,2));
  R2=V'*Z;
  Z=Z-V*R2;
  R1=R1+R2;
end

znrm=bvnorm(Z,2);

% orthogonalize Z
[Q,R2]=qr(bvscal(Z,1./znrm),0);

R2=bvscal(R2,znrm);

R(1:m,1:k)=R1;
R(m+1:m+k,1:k)=R2;

end
