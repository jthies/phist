n=5000;
cmplx=0;
m=100;
k=20;

relax=true;

% generate an orthogonal basis of m vectors
W1=randn(n,m) + i*cmplx*randn(n,m);
[V1,R1] = orthog(W1);

% k random new vectors W
Wj=randn(n,k)+i*cmplx*randn(n,k);

if (relax)
  disp('using classical GS');
else
  disp('using repeated classical GS');
end

% expand basis
[Vj,Rj]=orthog(V1,Wj,relax);

%test
disp('the following quantities should be SMALL');
disp('============================================');
disp('orthonormality of V after orthog() function:');
disp(norm(V1'*V1-eye(m)))
disp('correct QR decomp by orthog() function:');
disp(norm(V1*R1-W1))
%disp('orthogonality of new vectors Vj after orthog:');
%disp(norm(Vj'*Vj-eye(k)))

V=[V1 Vj];
R=[[R1;zeros(k,m)],Rj];

disp('orthogonality of expanded basis:');
disp(norm(V'*V-eye(m+k)))
disp('correctness of complete factorization:');
disp(norm(V*R-[W1 Wj]))

% test stability of orthogonalization on some hard test problem

aeps=1e-3;

V = [1 1 1;
     aeps aeps 0;
     aeps 0 aeps];

V(:,1)=V(:,1)./norm(V(:,1),2);

V1=V;
V2=V;
n=size(V,1);
k=size(V,2);
  
for i=2:k
  V1(:,i)=orthog(V1(:,1:i-1),V1(:,i),true);
  V2(:,i)=orthog(V2(:,1:i-1),V2(:,i),false);
end     

disp('"relaxed" orthog:');
disp(norm(V1'*V1-eye(k)));
disp('default orthog:');
disp(norm(V2'*V2-eye(k)))

% Householder QR, W=V*R
%{
%W=randn(m,k) + i*cmplx*randn(m,k);
%for j=1:k
%  [alpha,beta]=house(W(j:end,j));
%  v=W(j:end,j);
%  v(1)=v(1)-alpha;
%  W(j:end,j:end)=W(j:end,j:end) + beta*v*(v'*W(j:end,j:end))
%end
%}
