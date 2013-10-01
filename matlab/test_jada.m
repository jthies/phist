setpath;
%for debugging - make results reproducible
rand('seed',77);
randn('seed',42);
complex=0;

nx=8;
A=make_testmat(nx);
n=size(A,1);

k=1; % block size

% complex case
if (complex)
  A=A+0.25*sprandn(A)*i*complex;
end

% symmetric problem
%A=A+A'+2*speye(n);
%A2=A2+A2'+2*speye(n);

opts.verbose=true;
opts.debug=false;

% JD options
opts.arnoldi=false; % start with Arnoldi?
opts.switchTol=1.0e-3;
opts.numEigs=10;
opts.maxIter=250;
opts.tol=1.0e-6;
opts.minSpace=10;
opts.maxSpace=20;

opts.iterFun=@bgmres;


% options for the outer GMRES loop
opts.lsOpts.tol=1.0;
opts.lsOpts.maxIter=25;
opts.lsOpts.m=50; % restart/truncation parameter:

% symmetric
%s=1./sqrt(sum(abs(A),2));
%S=spdiags(s,0,n,n);
%A=S*A*S;

%s=1./sum(abs(A),2);
%S=spdiags(s,0,n,n);
%A=S*A; % row equilibration

%opts.precOp.compute=@comp_iluprec;


v0=randn(n,k);
if (complex)
  v0=v0+complex*i*randn(n,k);
end
v0=v0./norm(v0,2);

disp('JDQR for exterior eig(s).');
opts.target='LM';
[D,V,Q,R]=jdqre(A,v0,opts);

