function [lambda,v]=basic_jada(A,v0,opts,varargin)
%function [D,v]=basic_jada(A,v0,opts,varargin)
% basic Jacobi-Davidson for a single Eigenvalue
% near a shift (target). Works for non-Hermitian A.
% Ritz-Galerkin extraction is used, which makes the
% method less suitable for interior eigenvalues.
n=size(A,1);

global printOpts;

% some settings
maxIter=getopt(opts,'maxIter',30);
tol=getopt(opts,'tol',1.0e-6);

lsOpts.maxIter=10;
lsOpts.tol=0.1;

lsOpts = getopt(opts,'lsOpts',lsOpts);

target=getopt(opts,'target',0);
switchTol=getopt(opts,'switchTol',1.0e-3);

lsFun=getopt(opts,'iterFun',@bfgmres);
precOp=getopt(opts,'precOp',comp_idprec(A));


tau=target;
if (target==inf)
  tau=0;
end
% compute preconditioner and keep it fixed during iteration
precOp=precOp.compute(A-tau*speye(n),precOp);

% start Jacobi-Davidson
disp(sprintf('Jacobi-Davidon\n%s\t%s\t%s','iter','approx','resid'));

V=[];
t=v0; % starting guess
for m=1:maxIter
  t=mgs(V(:,1:m-1),t);
  V(:,m)=bvscal(t,1./bvnorm(t,2));
  AV(:,m)=A*V(:,m);
  %Hermitian case
%  H(1:m,m)=V(:,1:m)'*AV(:,m);
%  H(m,1:m-1)=H(1:m-1,m)';
   % non-Hermitian A
   H(1:m-1,m)=V(:,1:m-1)'*AV(:,m);
   H(m,1:m-1)=V(:,m)'*AV(:,1:m-1);
   H(m,m)=V(:,m)'*AV(:,m);
  % compute the largest eigenpair of H
  [S,D]=eig(H);
  [theta,i]=min(abs(diag(D)-target));
  theta=D(i,i);
  s=S(:,i);
  u=V*s;
  Au=AV*s;
  r=Au-theta*u;
  nrm=norm(r,2);
  disp(sprintf('%d\t%8.4e\t%8.4e',m,theta,nrm));
  
  if (nrm<tol)
    disp('converged!');
    lambda=theta;
    v=u;
    break;
  elseif (nrm<switchTol)
    disp('using RQI');
    shift=theta;
  else
    disp('using SI')
    shift=tau;
  end
  % solve approximately 
  % (I-uu')(A-theta*I)(I-uu')*t=-r
  % to get t \orth u

    % nested GMRES iterations
    op=comp_jada_op(A,shift,speye(n),u);

    disp('solve correction equation');
    disp(['preconditioner: ',precOp.label]);
    
    if (isfield(printOpts,'indent'))
      printOpts.indent=printOpts.indent+1;
    else
      printOpts.indent=1;
    end

    t0=t-(u'*t)*u;
    t0=t0./norm(t0,2);
    [t,flag,relres,iter,resvec]=lsFun(op,-r,t0,lsOpts,precOp);

    printOpts.indent=printOpts.indent-1;

end

