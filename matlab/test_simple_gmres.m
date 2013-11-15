setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
nx=8;

debug=false;

k=1; % number of rhs
% restart/truncation parameter:
m=20;

global printOpts;
printOpts.debug=debug;

A=make_testmat(nx,true);
n=size(A,1);

A2=make_testmat(nx)+sprandn(A)*i;

opts.debug=debug;
opts.tol=1e-6;
opts.maxIter=100;
opts.m=m;


% compare with a long GMRES series with
% up to 100MB of Krylov vectors
mlong=min(n,round(1e8/(8*n*k)));

s=1./sum(abs(A),2);
A=spdiags(s,0,n,n)*A;

s2=1./sum(abs(A),2);
A2=spdiags(s2,0,n,n)*A2;

% preconditioning:
%M=comp_iluprec(A);
%M=comp_diagprec(A);
M=comp_idprec(A);

XEX=randn(n,k);
XEX=bvscal(XEX,1./bvnorm(XEX));
XEX2=randn(n,k)+randn(n,k)*i;
XEX2=bvscal(XEX,1./bvnorm(XEX));
% for testing block variants
%XEX = repmat(XEX(:,1),1,k);
B=A*XEX;
B2=A2*XEX2;

X0=zeros(n,k);

opts.m=m;
disp(['real (B)GMRES(',int2str(m),')']);
[Xb,flag,relres,iter,resvec]=bgmres(A,B,X0,opts,M);

disp('');
disp(['real GMRES(',int2str(m),')']);
[Xb,flag,relres,iter,resvec]=simple_gmres(A,B,X0,opts,M);
disp('');
% COMPLEX
% preconditioning:
M=M.compute(A2,M);
opts.m=m;
disp(['complex (B)GMRES(',int2str(m),')']);
[Xb2,flag2,relres2,iter2,resvec2]=bgmres(A2,B2,X0,opts,M);
disp('');
disp(['complex (B)GMRES(',int2str(m),')']);
[Xb2,flag2,relres2,iter2,resvec2]=simple_gmres(A2,B2,X0,opts,M);
