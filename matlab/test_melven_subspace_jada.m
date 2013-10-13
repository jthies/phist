setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
complex=0;

nx=8;
A=make_testmat(nx);
n=size(A,1);

k=6; % block size

% complex case
if (complex)
  A=A+0.25*sprandn(A)*i*complex;
end


v0=randn(n,k);
if (complex)
  v0=v0+complex*i*randn(n,k);
end
v0=v0./norm(v0,2);

disp('subspace JaDa iteration for exterior eig(s).');

maxOuterIter = 10;
maxInnerIter = floor(80/k);
res_eps = 1.e-8;
resnorm_history_total = [];
for iOuter = 1:maxOuterIter
    fprintf('Outer iteration %d:', iOuter)
    [r,q,resnorm,resnorm_history,m] = melven_subspace_jada(A,v0,maxInnerIter,res_eps);
    resnorm_history_total = [resnorm_history_total resnorm_history];
    v0 = q;
    if( m < maxInnerIter )
        break;
    end
end
fprintf('Outer iterations: %d, total iterations: %d',iOuter, (iOuter-1)*maxInnerIter+m-1);
fprintf('\nEigenv. (MATLAB): %s', num2str(eigs(A,k)','%8.4g'));
fprintf('\n      Difference: '); fprintf(' %6.2g',abs(sort(diag(r),'descend')-sort(eigs(A,k),'descend')));
fprintf('\n');

