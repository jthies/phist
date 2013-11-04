function [D,V,Q,R]=jdqre(A,v0,opts)
%function [D,v]=jdqre(A,v0,opts,varargin)          
% single-vector JDQR for a few exterior eigenpairs.     
% JDQR method for kmax exterior eigenvalues, with       
% restart and deflation. Implementation follows         
% Eigentemplates p. 222 (alg. 7.18), and Sleijpen's     
% JDQR package.                                         
%                                                       
% input:                                                
%                                                       
% A - sparse matrix             
% v0 - starting vector          
% opts - options struct         
% M    - preconditioner.        
%                               
% options                       
% -------                       
%                               
% opts.numEigs [6]              
% opts.maxIter [200]            
% opts.tol [1e-6]               
% opts.target ['SM']            
% opts.maxSpace [25]            
% opts.minSpace [10]             
% opts.arno [true] use Arnoldi- 
%       iterations to start up, 
%       like in Sleijpens JDQR. 
% opts.switchTol [1e-2] [not used right now]
% opts.iterFun [bfgmres]        
% opts.lsOpts(passed to lsFun)  
% opts.precOp [idprec]          
% opts.verbose                  
% opts.debug                    
%                               
n=size(A,1);

global printOpts;
global count_MVM;
global A_operator;

A_operator=A; % for testing Sleijpens GMRES
count_MVM=0;

EXPAND=true;
SOLVE=true;

% some settings
numEigs=getopt(opts,'numEigs',6);
maxIter=getopt(opts,'maxIter',200);
tol=getopt(opts,'tol',1.0e-6);
mmax=getopt(opts,'maxSpace',max(25,2*numEigs));
mmin=getopt(opts,'minSpace',max(10,numEigs+4));

arno=getopt(opts,'arnoldi',true);

verbose=getopt(opts,'verbose',true);
debug=getopt(opts,'debug',false);

lsOpts.maxIter=10;
lsOpts.tol=0.1;

lsOpts = getopt(opts,'lsOpts',lsOpts);

target=getopt(opts,'target','LM');

switchTol=inf;
if (~arno)
  switchTol=1e-3;
end
switchTol=getopt(opts,'switchTol',switchTol);

lsFun=getopt(opts,'iterFun',@bfgmres);
precOp=getopt(opts,'precOp',comp_idprec(A));

tau=target;
if ischar(target)
  tau=0;
end

% compute preconditioner and keep it fixed during iteration
precOp=precOp.compute(A-tau*speye(n),precOp);

k=0; % number of converged Eigenpairs
it=0; % total number of iterations
m=0; % number of iterations since last restart
t=v0; % starting guess
V=[];
Q=zeros(n,0);
R=[];

if (arno)
  if (verbose)
    disp(sprintf('starting with %d steps of Arnoldi',mmin));
  end
  % construct H(1:mmin+1,1:mmin) and V(:,1:mmin+1)
  % such that A*V(:,1:mmin) = V*H
  arno_opts.maxIter=mmin;
  [V,H]=arnoldi(A,v0,arno_opts);
  size(H)
  AV=V*H;
  V=V(:,1:mmin);
  %AV=A*V;
  %M=V'*AV;
  M=H(1:mmin,1:mmin);
  m=mmin;
  EXPAND=false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start Jacobi-Davidson                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('Jacobi-Davidson\n%s\t%s\t%s\t\t%s','iter','m','approx','resid'));

mm=0;

while (k<numEigs && it<=maxIter)
  if (EXPAND)
    m=m+1;
    mm=mm+1;
    it=it+1;
    t=orthog([Q,V(:,1:m-1)],t);
    V(:,m)=t; %bvscal(t,1./bvnorm(t,2));
    AV(:,m)=apply_op(V(:,m),A);
    % Galerkin - non-Hermitian A
    M(1:m-1,m)=V(:,1:m-1)'*AV(:,m);
    M(m,1:m-1)=V(:,m)'*AV(:,1:m-1);
    M(m,m)=V(:,m)'*AV(:,m);
  
    if (debug)
      disp('new vector');
      V(1:5,m)
    end
  else
    EXPAND=true;
  end
  
  % sorted Schur-decomposition, M=STS', with
  % the eigenvalues on the diagonal of T sorted
  % according to the target. We take Sleijpens 
  % implementation for this, it first computes 
  % the Schur form using schur(...) (LAPACK) and
  % then reorders it using Givens rotations.
  
  %TODO - only sort complete Schur form if restart
  %       is impending, cf. args 3 and 4 of the
  %       function
  [S,T]=SortSchur(M,target);%m==mmax,mmin);
  if (debug)
    disp('sorted Schur (upper k x k block)');
    T(1:min(m,numEigs),1:min(m,numEigs))
  end
  theta=T(1,1);
  s=S(:,1);
  u=V*s;
  Au=AV*s;
  r=Au-theta*u;

  atil = Q'*r;
  rtil = r-Q*atil;
  nrm=norm(rtil);

  print_eigs_iter(it,m,nrm,theta);

  % deflate converged eigenpairs
  %while (norm(r)<=tol)
  while (nrm<=tol)
    R=[R, atil;
       zeros(1,k), theta];
    Q=[Q,u];
    k=k+1;
    mm=0;
    if (verbose)
      disp(['eigenvalue ',int2str(k),' (',num2str(theta),') is converged.']);
      %more on;
    end %if
    if (debug)
      disp('updated Schur form (R):');
      R
    end
    if (k==numEigs)
      SOLVE=false;
      break;
    end

%    if (abs(imag(theta))>tol)
%      t=imag(u/sign(max(u)));
%      if (norm(t)>tol)
%        t=orthog(Qschur,t); 
%        EXPAND=(norm(t)>sqrt(tol)); 
%      end
%    end
    S=S(:,2:m);
    M=T(2:m,2:m);
    m=m-1;
    V=V*S;
    AV=AV*S;
    u=V(:,1);
    Au=AV(:,1);
    theta=M(1,1);
    r=Au-theta*u;
    atil=Q'*r;
    rtil=r-Q*atil;
    nrm=norm(rtil);

    [S,T]=SortSchur(M,target);%,m==mmax,mmin);

  end % while (deflate)

  % restart if necessary
  if m>=mmax
    if (verbose)
      disp('restart JDQR');
    end
    m=mmin;
    S=S(:,1:m);
    M=T(1:m,1:m);
    V=V*S;
    AV=AV*S;
    V=V(:,1:mmin);
    AV=AV(:,1:mmin);
    SOLVE=true;
    EXPAND=true;
  end
  if (SOLVE)
    % maintain orthogonality against
    % all converged vectors (Q) and the
    % new one:
    Qtil=[Q,u];
    if (nrm<switchTol)
      if (verbose)
        disp('using RQI');
      end
      shift=theta;
    else
      if (verbose)
        disp('using SI')
      end
      shift=tau;
    end
    % solve approximately 
    % (I-uu')(A-theta*I)(I-uu')*t=-r
    % to get t \orth u (u ^= Qtil here)
    if (strcmp(lsFun,'direct'))
      A_aug = [A-shift*speye(n), Qtil;
                Qtil'        , zeros(size(Qtil,2))];
      b_aug = [rtil;zeros(size(Qtil,2),size(rtil,2))];
      t_aug= A_aug\b_aug;
      t=t_aug(1:size(Qtil,1),:);
    else
      op=comp_jada_op(A,shift,speye(n),Qtil);
      if (isfield(printOpts,'indent'))
        printOpts.indent=printOpts.indent+1;
      else
        printOpts.indent=1;
      end
      % TODO - is this necessary?
      precOp=precOp.compute(A-shift*speye(n),precOp);
      t0=zeros(size(rtil));
      lsOpts.tol=max(tol,1/2.^max(1,mm));
      if (verbose)
        disp(['inner conv tol: ',num2str(lsOpts.tol)]);
      end
      
      [t,flag,relres,iter,resvec] = ...
        lsFun(op,rtil,t0,lsOpts,precOp);

      printOpts.indent=printOpts.indent-1;
    end % direct or user-supplied
    EXPAND=true;    
    if (verbose)
      disp(['(Qtil,t)=',num2str(norm(Qtil'*t))]);
      disp(['(A-sI)t+rtil=',num2str(norm(A*t-shift*t+rtil))]);         
    end
  else
    EXPAND=false;
    SOLVE=true;
  end %SOLVE
end % while

[z,D]=Jordan(R);
V=Q*z;

if (verbose)
  str=sprintf('number of iterations: %d',it);
  if (arno)
    str=[str, sprintf('(+%d Arnoldi steps)',mmin)];
  end
  disp(str);
  disp('computed Eigenvalues:');
  Resid=bvnorm(A*V-V*D);
  for i=1:k
    print_eigs_iter(i,[],Resid(i),D(i,i));
  end
if (verbose)
  disp(sprintf('Number of MVMs: %d',count_MVM));
end  

end
