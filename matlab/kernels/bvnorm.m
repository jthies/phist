function nrm=bvnorm(V,varargin)
% function nrm=bvnorm(V,type)
% function nrm=bvnorm(V) assumes type=2
% computes norm(V(:,j),type) for all columns j of V

if (nargin==1)
  type=2;
else
  type=varargin{1};
end

k=size(V,2);
nrm=zeros(1,k);
for i=1:k
  nrm(i) = norm(V(:,i),type);
end

end
