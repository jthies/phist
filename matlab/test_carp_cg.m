clear;
setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
nx=8;

debug=false;

k=1; % number of rhs
% restart param for GMRES
m=25;
tol=1.0e-12;
maxRestarts=20;
maxIt=500;

%A=make_testmat(nx,false);
A=mmread('gordon/matrix.mtx');
n=size(A,1);

A=-spones(A);
A=A+speye(n)*7;

% this matrix is similar to our Graphene matrices:
%A=-spones(A);
%d=randn(n,1);
%d=d./(2*max(abs(d)));
%A=spdiags(d,0,A);

sigma=0; % shift

xex=load('gordon/xex.txt');
b=load('gordon/rhs.txt');

x0=zeros(n,1);

% smallest eigenpairs for playing with explicit deflations
nd=4;
%[Vr,Lambda1]=eigs(A,nd,'SM');
%[Vl,Lambda2]=eigs(A',nd,'SM');
[V,D]=eigs(A,nd,'SM');
nd=size(D,1);

[U,R]=qr(V,0);
T=(R*D)/R;

disp('=============================');
disp('unpreconditioned GMRES');
disp('=============================');


[x1,flag1,relres1,iter1,resvec1] = gmres(A, b, m, tol, maxRestarts, [], [], x0);

disp(sprintf('total iters %d, restarts %d, relres %e\n',length(resvec1), iter1(1)-1, relres1));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x1-sigma*x1-b),norm(x1-xex)));

if (flag1~=0)
  disp(sprintf('non-zero return flag=%d',flag1));
end

disp('=============================');
disp(sprintf('GMRES orth. %d modes',nd));
disp('=============================');

[x2,flag2,relres2,iter2,resvec2] = gmres(A, b, m, tol, ...
        maxRestarts, [], @apply_deflprec, x0,speye(n),T,U);

disp(sprintf('total iters %d, restarts %d, relres %e\n',length(resvec2), iter2(1)-1, relres2));
%disp(sprintf('expl. resid: %e, error %e\n',norm(A*x2-sigma*x2-b),norm(x2-xex)));

if (flag2~=0)
  disp(sprintf('non-zero return flag=%d',flag2));
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

opts.omega=1.5;
[x4,flag4,relres4,iter4,resvec4] = carp_cg(A, b, x0, opts);
disp(sprintf('CARP-CG (omg=%3.1f): iters %d, relres %e\n',opts.omega,iter4, relres4));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x4-sigma*x4-b),norm(x4-xex)));

if (flag4~=0)
  disp(sprintf('non-zero return flag=%d',flag4));
end

disp(sprintf('CARP-CG with deflation of %d vectors',nd));

opts.Precond=comp_deflprec(speye(n),T,U);

opts.omega=1.0;
opts.tol=tol;
opts.maxIter=maxIt;

[x5,flag5,relres5,iter5,resvec5] = carp_cg(A, b, x0, opts);
disp(sprintf('CARP-CG: iters %d, relres %e\n',iter5, relres5));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x5-sigma*x5-b),norm(x5-xex)));

if (flag5~=0)
  disp(sprintf('non-zero return flag=%d',flag5));
end

opts.omega=1.5;
[x6,flag6,relres6,iter6,resvec6] = carp_cg(A, b, x0, opts);
disp(sprintf('CARP-CG (omg=%3.1f): iters %d, relres %e\n',opts.omega,iter6, relres6));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x6-sigma*x6-b),norm(x6-xex)));

if (flag6~=0)
  disp(sprintf('non-zero return flag=%d',flag6));
end

