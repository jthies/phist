A=make_testmat(16);
n=size(A,1);
k=4;
V=randn(n,4);
V=orthog(V);
sigma=0.3; % shift
I=speye(n);
op=comp_jada_op(A,sigma,I,V);

xex=randn(n,1);
x0=randn(n,1);

x0=orthog(V,x0);
xex=orthog(V,xex);

b=apply_op(xex,op);

opts.m=999;
opts.tol=1e-6;
opts.maxIter=999;

M=comp_iluprec(A-sigma*I);
Morth=comp_jada_op(M,0,I,V);
[x,flag,relres,iter,resvec]=bgmres(op,b,x0,opts,Morth);
disp('number of iterations:');
disp(iter);

disp('rel. error norm:');
disp(norm(x-xex)./norm(b));

disp('orthogonality wrt. V:');
disp(norm(V'*x));
