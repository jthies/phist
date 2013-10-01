function [X,Jordan]=Jordan(S)
% [X,J]=JORDAN(S)
%   For S  k by k upper triangular matrix with ordered diagonal elements,
%   JORDAN computes the Jordan decomposition.
%   X is a k by k matrix of vectors spanning invariant spaces.
%   If J(i,i)=J(i+1,i+1) then X(:,i)'*X(:,i+1)=0.
%   J is a k by k matrix, J is Jordan such that S*X=X*J.
%   diag(J)=diag(S) are the eigenvalues

% coded by Gerard Sleijpen, Januari 14, 1998

k=size(S,1); X=zeros(k);
if k==0, Jordan=[]; return, end

%%% accepted separation between eigenvalues:
delta=2*sqrt(eps)*norm(S,inf); delta=max(delta,10*eps);

T=eye(k); s=diag(S); Jordan=diag(s);
for i=1:k
  I=[1:i]; e=zeros(i,1); e(i,1)=1;
  C=S(I,I)-s(i,1)*T(I,I); C(i,i)=1;
  j=i-1; q=[]; jj=0;
  while j>0
    if abs(C(j,j))<delta, jj=jj+1; j=j-1; else, j=0; end
  end
  q=X(I,i-jj:i-1);
  C=[C,T(I,I)*q;q',zeros(jj)];
  q=C\[e;zeros(jj,1)]; nrm=norm(q(I,1));
  Jordan(i-jj:i-1,i)=-q(i+1:i+jj,1)/nrm;
  X(I,i)=q(I,1)/nrm;
end

return

