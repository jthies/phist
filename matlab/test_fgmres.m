%for debugging - make results reproducible
%rand("seed",42);
%randn("seed",42);

nx=16;

A=make_testmat(nx);
A2=make_testmat(nx)+ 0.25*i*sprandn(A);
n=size(A,1);
k=2; % number of rhs
m=4; % restart/truncation parameter

opts.outer.debug=false;
opts.inner.debug=false;

% compare with a long GMRES series with
% up to 1GB of Krylov vectors
mlong=round(min(n,1e9/(8*n*k)));

% options for the outer GMRES loop
opts.outer.tol=1.0e-6;
opts.outer.maxIter=100;
opts.outer.m=m; % restart/truncation parameter:

% nested GMRES iterations
opts.inner.tol=0.1;
opts.inner.maxIter=5;
opts.inner.m=5; % restart/truncation parameter:

s=1./sum(abs(A),2);
A=spdiags(s,0,n,n)*A;

s2=1./sum(abs(A),2);
A2=spdiags(s2,0,n,n)*A2;

XEX=randn(n,k);
XEX=bvscal(XEX,1./bvnorm(XEX));
% for testing block variants
%XEX = repmat(XEX(:,1),1,k);
B=A*XEX;
B2=A2*XEX;

X0=zeros(n,k);

M.opts=opts.inner;
M=comp_iterprec(A,M);

disp(['real (B)FGMRES(',int2str(mlong),')']);
opts.outer.m=mlong;
[Xa,flag,relres,iter,resvec]=bfgmres(A,B,X0,opts.outer,M);
disp(['real (B)FGMRES(',int2str(m),')']);
opts.outer.m=m;
[Xb,flag,relres,iter,resvec]=bfgmres(A,B,X0,opts.outer,M);
disp(['real (B)DQ-GMRES(',int2str(m),')']);
[Xc,flag,relres,iter,resvec]=bdqgmres(A,B,X0,opts.outer,M);

disp(['complex (B)FGMRES(',int2str(mlong),')']);
[Xa2,flag2,relres2,iter2,resvec2]=bfgmres(A2,B2,X0,opts.outer,M);
disp(['complex (B)FGMRES(',int2str(m),')']);
[Xb2,flag2,relres2,iter2,resvec2]=bfgmres(A2,B2,X0,opts.outer,M);
disp(['complex (B)DQ-GMRES(',int2str(m),')']);
[Xc2,flag2,relres2,iter2,resvec2]=bdqgmres(A2,B2,X0,opts.outer,M);

% GMRES reference implementation (with left preconditioning and MGS)
%[X,err,iter,flag]=template_gmres(A,X0,B,speye(n,n),m,maxIt,tol);
%[X2,err2,iter2,flag2]=template_gmres(A2,X0,B2,speye(n,n),m,maxIt,tol);
