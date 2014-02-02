setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
nx=8;

debug=false;

k=1; % number of rhs
% restart param for GMRES
m=25;
tol=1.0e-4;
maxRestarts=20;
maxIt=500;

A=make_testmat(nx,true);
n=size(A,1);

% this matrix is similar to our Graphene matrices:
A=-spones(A);
d=randn(n,1);
d=d./(2*max(abs(d)));
A=spdiags(d,0,A);

sigma=0; % shift

xex=ones(n,1);
xex=xex/norm(xex);
b=A*xex-sigma*xex;
x0=zeros(n,1);

disp('=============================');
disp('diagonal preconditioned GMRES');
disp('=============================');
%B=speye(n);
%M=spdiags(full(sum(abs(A-sigma*B),2)),0,n,n);
%A=A*inv(M);
%op=comp_shifted_op(M\A,-sigma,M\B);
%[x1,flag1,relres1,iter1,resvec1] = gmres(@apply_op, M*b, m, tol, maxIt, [], [], x0,op);
[x1,flag1,relres1,iter1,resvec1] = gmres(A, b, m, tol, maxRestarts, [], [], x0);

disp(sprintf('total iters %d, restarts %d, relres %e\n',length(resvec1), iter1(1)-1, relres1));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x1-sigma*x1-b),norm(x1-xex)));

if (flag1~=0)
  disp(sprintf('non-zero return flag=%d',flag1));
end

%disp('=====================');
%disp('unpreconditioned CGNR');
%disp('=====================');

%op=comp_AAt_op(A,sigma);

%[y2,flag2,relres2,iter2,resvec2] = pcg(@apply_op, b, tol, maxIt, [], [], x0,op);
%x2=A'*y2;
%disp(sprintf('CGNR: iters %d, relres %e\n',iter2, relres2));
%disp(sprintf('expl. resid: %e, error %e\n',norm(A*x2-sigma*x2-b),norm(x2-xex)));

%if (flag2~=0)
%  disp(sprintf('non-zero return flag=%d',flag2));
%end


opts.omega=1.0;
opts.tol=tol;
opts.maxIter=maxIt;
[x3,flag3,relres3,iter3,resvec3] = carp_cg(A, b, x0, opts);
disp(sprintf('CARP-CG: iters %d, relres %e\n',iter3, relres3));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x3-sigma*x3-b),norm(x3-xex)));

if (flag3~=0)
  disp(sprintf('non-zero return flag=%d',flag3));
end

opts.omega=1.9;
[x4,flag4,relres4,iter4,resvec4] = carp_cg(A, b, x0, opts);
disp(sprintf('CARP-CG (omg=%3.1f): iters %d, relres %e\n',opts.omega,iter4, relres4));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x4-sigma*x4-b),norm(x4-xex)));

if (flag4~=0)
  disp(sprintf('non-zero return flag=%d',flag4));
end

