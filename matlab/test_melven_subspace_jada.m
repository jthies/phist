setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
complex=1;

nx=8;
A=make_testmat(nx);
n=size(A,1);

k=8; % block size

% complex case
if (complex)
  A=A+0.25*sprandn(A)*i*complex;
end

generalized = false;
if( generalized )
  % create spd. B
  B = sprandsym(n, 0.001,0.1,2);
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
maxIter = 400;
res_eps = 1.e-8;
resnorm_history_total = [];
[r,q,resnorm,resnorm_history,m,restarts] = melven_subspace_jada(A,B,v0,maxIter,minBas,maxBas,res_eps);
fprintf('total iterations: %d, restarts: %d',m,restarts);
fprintf('\nEigenv. (MATLAB): %s', num2str(eigs(A,B,k)','%8.4g'));
fprintf('\n      Difference: '); fprintf(' %6.2g',abs(sort(diag(r),'descend')-sort(eigs(A,B,k),'descend')));
fprintf('\n');

