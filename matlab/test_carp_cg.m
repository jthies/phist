clear all;
close all;
setpath;

mex -largeArrayDims krylov/dkswp2.c
mex -largeArrayDims krylov/nrms_ai2.c

%for debugging - make results reproducible
%rand('seed',77);
%randn('seed',42);

debug=false;

k=1; % number of rhs

% for GMRES:
m=50; % max size of Krylov basis 
maxRestarts=200;
% for all methods
maxIt=10000;
tol=1.0e-12;

mpath='/hpc_data/essex/WPT/';
mpath2='/home/thie_jo/work/fredwubs/hymls/testSuite/data/Graphene/';

matrix_fmt='mm';
matrices={[mpath2,'128x64/jac.mtx'],
          [mpath2,'256x128/jac.mtx'],
          [mpath2,'512x256/jac.mtx']};

%matrix_fmt='mat';
%matrices={[mpath,'graphen/graphen21x4000.mat'],
%          [mpath,'graphen/graphen22x8000.mat'],
%          [mpath,'graphen/graphen21x40000.mat']};
load([mpath,'shifts/graphene_shifts.mat']);
shifts=graphene_shifts;


%matrix_fmt='mat';
%matrices={[mpath,'lap_cit/LAP_CIT_396.mat'],
%          [mpath,'lap_cit/LAP_CIT_1059.mat'],
%          [mpath,'lap_cit/LAP_CIT_3084.mat'],
%          [mpath,'lap_cit/LAP_CIT_4470.mat'],
%          [mpath,'lap_cit/LAP_CIT_6752.mat'],
%          [mpath,'lap_cit/LAP_CIT_8843.mat']}

skipGMRES=true;
for matrixID=1:length(matrices)

if strcmp(matrix_fmt,'mat')
  tmp=load(matrices{matrixID});
  A=tmp.S;
elseif strcmp(matrix_fmt,'mm')
  A=mmread(matrices{matrixID});
else
  error('matrix format "',matrix_fmt,'" not known');
end
n=size(A,1);
I=speye(n);
xex=randn(n,1)./sqrt(n);
x0=zeros(n,1);

% generate 4-color ordering for Graphene (special case)
nx=128*2^(matrixID-1);
ny=nx/2;
ord=color_graphene(nx,ny);
A=A(ord,ord);
xex=xex(ord);

for shiftID=1:length(shifts)
  sigma=shifts(shiftID);
  disp(['TEST CASE: ',matrices{matrixID},' sigma=',num2str(sigma)]);
  opts=struct;

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
[x3,flag3,relres3,iter3,resvec3,Tlan3] = carp_cg(A, b, x0, opts);
toc
disp(sprintf('CARP-CG: iters %d, relres %e\n',iter3, relres3));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x3-sigma*x3-b),norm(x3-xex)));

if (flag3~=0)
  disp(sprintf('non-zero return flag=%d',flag3));
end

% try some overrelaxation
%{
omg2=0.85;
opts.omega=omg2;
tic;
[x4,flag4,relres4,iter4,resvec4,Tlan4] = carp_cg(A, b, x0, opts);
toc
disp(sprintf('CARP-CG (omg=%3.1f): iters %d, relres %e\n',opts.omega,iter4, relres4));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x4-sigma*x4-b),norm(x4-xex)));

if (flag4~=0)
  disp(sprintf('non-zero return flag=%d',flag4));
end
%}

% try deflation
%{
nseed=16;
I=speye(n);
%d=nrms_ai2(A-sigma*I);
%D=spdiags(d,0,n,n);
%[V,Lambda]=eigs((A-real(sigma)*I)/D,nseed,'SM');
%[V,Lambda]=eigs(A,nseed,real(sigma));
V=orth(randn(n,nseed));
opts.omega=1.0;
%z=zeros(nseed,1);
%opts.massmat=[I,sparse(n,nseed);sparse(nseed,n+nseed)];
opts.seedspace=V;
%opts.Precond=comp_deflprec(I,V'*(A*V)-sigma*eye(nseed),V);
tic;
[x5,flag5,relres5,iter5,resvec5,Tlan5] = carp_cg(A, b, x0, opts);
%[y5,flag5,relres5,iter5,resvec5] = carp_cg([A, V; V',sparse(nseed,nseed)], ...
%        [b-V*(V'*b);z], [x0;z],opts);
%x5=y5(1:n);        
toc
disp(sprintf('CARP-CG (defl %d): iters %d, relres %e\n',nseed,iter5, relres5));
disp(sprintf('expl. resid: %e, error %e\n',norm(A*x5-sigma*x5-b),norm(x5-xex)));

if (flag5~=0)
  disp(sprintf('non-zero return flag=%d',flag5));
end
%}

figure(matrixID);
hleg=legend;
if isempty(hleg)
  kleg=0;
else
  leg_ud=get(hleg,'UserData');
  leg=leg_ud.lstrings;
  kleg=size(leg,1);
end
if exist('x1')
  semilogy(resvec1,'k-.');
  kleg=kleg+1;
  leg{kleg}=sprintf('GMRES(%d), matrix %d, shift %d',m,matrixID,shiftID);
  hold all;
end
if exist('x2')
  semilogy(resvec3,'k-+');
  kleg=kleg+1;
  leg{kleg}=sprintf('GMRES(%d) + ILUTP, matrix %d, shift %d',m,matrixID, shiftID);
  if ~ishold
    hold all;
  end
end
if exist('x3')
  semilogy(resvec3,'k-');
  kleg=kleg+1;
  leg{kleg}=sprintf('CARP-CG, matrix %d, shift %d',matrixID, shiftID);
  if ~ishold
    hold all;
  end
end
if exist('x4')
  semilogy(resvec4,'k--');
  kleg=kleg+1;
  leg{kleg}=sprintf('CARP-CG(%f), matrix %d, shift %d',omg2, matrixID, shiftID);
  if ~ishold
    hold all;
  end
end
if exist('x5')
  semilogy(resvec5,'k:');
  kleg=kleg+1;
  leg{kleg}=sprintf('CARP-CG(defl %d), matrix %d, shift %d',nseed, matrixID, shiftID);
  if ~ishold
    hold all;
  end
end
legend(leg);
%title(matrices{matrixID});
%hold off;
end


%ritz = sort(eig(full(Tlan3)));

end % matrixID (test cases)