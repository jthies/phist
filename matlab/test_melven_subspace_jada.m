setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
complex=1;

nx=8;
A=make_testmat(nx);
n=size(A,1);


k = 4;
nEig = 20;
res_eps = 1.e-8;
generalized = false;

%% for generalized
%k=4; % block size
%nEig = 10;
%res_eps = 1.e-3;
%generalized = true;
%%

% complex case
if (complex)
  A=A+0.25*sprandn(A)*1i*complex;
end

if( generalized )
  % create spd. B
  B = sprandsym(n, 0.001,0.5,2);
  % print cond of B
  fprintf('cond(B) = %f\n', condest(B));
else
  B = speye(n);
end

v0=randn(n,k);
if (complex)
  v0=v0+complex*1i*randn(n,k);
end
v0=v0./norm(v0,2);

disp('subspace JaDa iteration for exterior eig(s).');

minBas = 2*k;
maxBas = 80;
arnoldiIter = 2;
maxIter = 400;
resnorm_history_total = [];
[r,q,resnorm,lambda_history,resnorm_history,m,restarts] = melven_subspace_jada(A,B,v0,nEig,arnoldiIter,maxIter,minBas,maxBas,res_eps);
fprintf('total iterations: %d, restarts: %d',m,restarts);
fprintf('\nEigenv. (MATLAB): %s', num2str(eigs(A,B,nEig)','%8.4g'));
fprintf('\n      Difference: '); fprintf(' %6.2g',abs(sort(diag(r),'descend')-sort(eigs(A,B,nEig),'descend')));
fprintf('\n');

