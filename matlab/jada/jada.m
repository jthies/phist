function [lambda,v]=jada(A,v0,opts,varargin)
%                                                               
%function [D,v]=jada(A,v0,opts,varargin)                        
% basic Jacobi-Davidson for opts.nev exterior Eigenpairs        
% or Eigenpairs near a target. Works for non-Hermitian A.       
% Ritz-Galerkin extraction is used for exterior eigs, harmonic  
% Ritz for interior ones.                                       
%                                                               
% options                                                       
% ~~~~~~~                                                       
%                                                               
% opts.numEigs
% opts.target [inf] a real or complex number. 'inf' means exte- 
%               rior eigenvalues are sought.                    
% opts.maxIter                                                  
% opts.tol                                                      
% opts.switchTol - at which residual norm to switch to RQI      
% opts.iterFun [@bgmres] - which iterative solver to use        
% opts.precOp [no preconditioning] - preconditioning operator   
%               for the iterative linear solver.                
% opts.outer - options to pass to the iterative linear solver   
%       opts.outer.tol: if set to 1 (the default), the value    
%       2^-j is used in Jada-iteration j.
% opts.inner - can be used for nested iterative solvers         
n=size(A,1);
v=v0;
k=size(v,2); % block size (must be 1 up to now)
global printOpts;

precComputed=false;

% some default settings
nev=1; % number of eigenpairs desired
maxIter=30;
tol=opts.tol;
outerGmresOpts.maxIter=25;
outerGmresOpts.tol=1.0;
innerGmresOpts.maxIter=10;
innerGmresOpts.tol=0.1;
target=inf;
switchTol=1.0e-3;

gmresFun=@bfgmres;
precOp=comp_idprec(A);

if (isfield(opts,'numEigs'))
  nev = opts.numEigs;
end
if (isfield(opts,'maxIter'))
  maxIter = opts.maxIter;
end
if (isfield(opts,'outer'))
  outerGmresOpts = opts.outer;
end
if (isfield(opts,'inner'))
  innerGmresOpts = opts.inner;
end
if (isfield(opts,'iterFun'))
  gmresFun=opts.iterFun;
end
if (isfield(opts,'precOp'))
  precOp=opts.precOp;
end
% if target is not specified, use theta from the start
% (gives largest magnitude eigenvalue)
if (isfield(opts,'target'))
  target=opts.target;
end
if (isfield(opts,'switchTol'))
  switchTol=opts.switchTol;
end

%if (target==inf)
%  switchTol=inf; % use RQI directly
%end
adaptTol=1.0;
if (outerGmresOpts.tol==1.0)
  adaptTol=0.5;
end
% start Jacobi-Davidson
disp(sprintf('Jacobi-Davidon\n%s\t%s\t%s','iter','approx','resid'));

V=[];
t=v0; % starting guess
for m=1:maxIter
  t=mgs(V(:,1:m-1),t);
  V(:,m)=bvscal(t,1./bvnorm(t,2));
  W(:,m)=A*V(:,m);
  % Hermitian case
  %H(1:m,m)=V(:,1:m)'*W(:,m);
  %H(m,1:m-1)=H(1:m-1,m)';
  % non-Hermitian A
  H(1:m-1,m)=V(:,1:m-1)'*W(:,m);
  H(m,1:m-1)=V(:,m)'*W(:,1:m-1);
  H(m,m)=V(:,m)'*W(:,m);
  % compute the largest eigenpair of H
  % (Ritz-Galerkin extraction)
  if (target==inf)
    [S,D]=eig(H);
    [theta,i]=max(abs(diag(D)));
    theta=D(i,i);
  else
    %disp('warning - using Ritz-Galerkin for interior eigs');
    [S,D]=eig(H);
    [theta,i]=min(abs(diag(D)-target));
    theta=D(i,i);
  end
  s=S(:,i);
  u=V*s;
  Au=W*s;
  r=Au-theta*u;
  nrm=norm(r,2);
  disp(sprintf('%d\t%8.4f\t%8.4f',m,theta,nrm));
  
  if (nrm<tol)
    lambda=theta;
    v=u;
    break;
  elseif (nrm<switchTol||target==inf)
    shift=theta;
  else
    shift=target;
  end
  % solve approximately 
  % (I-uu')(A-theta*I)(I-uu')*t=-r
  % to get t \orth u
  
  % compute preconditioner and keep it fixed during iteration
  if (~precComputed)
    precOp=precOp.compute(A-shift*speye(n),precOp);
    precComputed=true;
  end


    % nested GMRES iterations
    op=comp_jada_op(A,shift,speye(n),u);
    
    if (isfield(printOpts,'indent'))
      printOpts.indent=printOpts.indent+1;
    else
      printOpts.indent=1;
    end

    t0=t-(u'*t)*u;
    t0=t0./norm(t0,2);
    outerGmresOpts.tol = outerGmresOpts.tol*adaptTol;
    [t,flag,relres,iter,resvec]=gmresFun(op,-r,t0,outerGmresOpts,...
          precOp);

    printOpts.indent=printOpts.indent-1;

end

