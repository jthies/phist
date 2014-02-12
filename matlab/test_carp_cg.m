clear all;
close all;
setpath;

mex -largeArrayDims krylov/dkswp2.c
mex -largeArrayDims krylov/nrms_ai2.c

%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);

debug=false;

k=1; % number of rhs

% for GMRES:
m=50; % max size of Krylov basis 
maxRestarts=200;
% for all methods
maxIt=10000;
tol=1.0e-12;

mpath='/hpc_data/essex/WPT/';
matrices={[mpath,'graphen/graphen21x4000.mat'],
          [mpath,'graphen/graphen22x8000.mat'],
          [mpath,'graphen/graphen21x40000.mat']};

%matrices={[mpath,'lap_cit/LAP_CIT_396.mat'],
%          [mpath,'lap_cit/LAP_CIT_1059.mat'],
%          [mpath,'lap_cit/LAP_CIT_3084.mat'],
%          [mpath,'lap_cit/LAP_CIT_4470.mat'],
%          [mpath,'lap_cit/LAP_CIT_6752.mat'],
%          [mpath,'lap_cit/LAP_CIT_8843.mat']}


% solve (A-sigma*I)x=b
sigma=0.008+1.0e-4i;

skipGMRES=true;

for matrixID=1:length(matrices)

disp(['TEST CASE: ',matrices{matrixID}]);

tmp=load(matrices{matrixID});
A=tmp.S;
n=size(A,1);
I=speye(n);
xex=randn(n,1)./sqrt(n);
x0=zeros(n,1);


b=A*xex-sigma*xex;
%b=b-V*V'*b;
%
if (skipGMRES==false)

disp('=============================');
disp(sprintf('unpreconditioned GMRES(%d)',m));
disp('=============================');
%
tic;
[x1,flag1,relres1,iter1,resvec1] = gmres(A-sigma*I, b, m, tol, maxRestarts, [], [], x0);
toc
disp(sprintf('total iters %d, restarts %d, relres %e\n',length(resvec1), iter1(1)-1, relres1));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x1-sigma*x1-b),norm(x1-xex)));

if (flag1~=0)
  disp(sprintf('non-zero return flag=%d',flag1));
end
%
%
disp('=============================');
disp(sprintf('ILU-preconditioned GMRES(%d)',m));
disp('=============================');
tic;
iluOpts.type='ilutp';
iluOpts.droptol=0.01;
iluOpts.milu='off';
iluOpts.udiag=1;
iluOpts.thresh=1; % [0..1], 1: no pivoting
[L,U]=ilu(A-sigma*I,iluOpts);
%[L,U]=ilu(A); % ILU(0)
disp('prec setup');
toc
tic;
[x2,flag2,relres2,iter2,resvec2] = gmres(A-sigma*I, b, m, tol, maxRestarts, L, U, x0);
toc

disp(sprintf('total iters %d, restarts %d, relres %e\n',length(resvec2), iter2(1)-1, relres2));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x2-sigma*x2-b),norm(x2-xex)));

if (flag2~=0)
  disp(sprintf('non-zero return flag=%d',flag2));
end

end % skipGMRES?
%
opts.sigma=sigma;
opts.omega=1.0;
opts.tol=tol;
opts.maxIter=maxIt;
tic;
[x3,flag3,relres3,iter3,resvec3] = carp_cg(A, b, x0, opts);
toc
disp(sprintf('CARP-CG: iters %d, relres %e\n',iter3, relres3));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x3-sigma*x3-b),norm(x3-xex)));

if (flag3~=0)
  disp(sprintf('non-zero return flag=%d',flag3));
end

% try some overrelaxation
opts.omega=1.33;
tic;
[x4,flag4,relres4,iter4,resvec4] = carp_cg(A, b, x0, opts);
toc
disp(sprintf('CARP-CG (omg=%3.1f): iters %d, relres %e\n',opts.omega,iter4, relres4));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x4-sigma*x4-b),norm(x4-xex)));

if (flag4~=0)
  disp(sprintf('non-zero return flag=%d',flag4));
end
%

kleg=0;
figure(matrixID);
if exist('x1')
  semilogy(resvec1,'k-');
  kleg=kleg+1;
  leg{kleg}=sprintf('GMRES(%d)',m);
  hold on;
end
if exist('x2')
  semilogy(resvec3,'k-.');
  kleg=kleg+1;
  leg{kleg}=sprintf('GMRES(%d) + ILUTP',m);
  hold on;
end
if exist('x3')
  semilogy(resvec4,'k:');
  kleg=kleg+1;
  leg{kleg}=sprintf('CARP-CG');
  hold off;
end
if exist('x4')
  semilogy(resvec4,'k:');
  kleg=kleg+1;
  leg{kleg}=sprintf('CARP-CG(%f)',opts.omega);
  hold off;
end

legend(leg);
title(matrices{matrixID});
end
