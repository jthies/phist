k=5;
m=10; % number of steps

A=make_testmat(4);
n=size(A,1);
m=n;
% generate k orthogonal starting vectors
v0=randn(n,k);

opts.maxIter = m;

[V,H]=arnoldi(A,v0,opts);
m=size(V,2)./k-1;

disp('as usual, these numbers should be small.');
disp('quality of factorization:');
disp(norm(A*V(:,1:m*k) - V*H));
disp('orthogonality of basis');
disp(norm(V'*V-eye((m+1)*k)));
