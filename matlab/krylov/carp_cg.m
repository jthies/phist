function [x,flag,relres,iter,resvec,Tlan]=carp_cg(A,b,x0,opts)
%                                                               
% [x,flag,relres,iter,resvec, Tlan]=carp_cg(A,b,x0,opts)        
%                                                               
% CGMN algorithm (Bjoerck & Elfving 1979). Returns the typical  
% output of iterative methods in matlab (cf. 'help pcg'), and   
% optionally computes the tridiagonal Lanczos matrix so that    
% eigenvalue information can be obtained a posteriori.          
%                                                        
% options                                                
% ~~~~~~~                                                
%                                                        
% opts.tol - convergence tolerance                       
% opts.maxIter - max number of iterations                
% opts.omega - relaxation parameter                      
% opts.Precond - preconditioning operator                
% opts.sigma - shift, solve (A-sigma*I)x=b instead.      
%              in contrast to A, sigma may be complex.   
% opts.seedSpace (experimental) - attempt to deflate a   
% number of approximate eigenvectors.                      
%
verbose=getopt(opts,'verbose',true);

itprint=1;
itcheck=1;

errormode = false;
correction_needed = false;


flag=0;
relres=1.0;
iter = 0; % total iterations
resvec=[];

wantT=(nargout>=6);
n=size(A,1);

cfreq = getopt(opts,'cfreq',100);          %correction frequency
tol=getopt(opts,'tol',1e-8);
maxIter=getopt(opts,'maxIter',300);
omega=getopt(opts,'omega',1.7);
M=getopt(opts,'Precond',speye(n));
seedSpace=getopt(opts,'seedspace',[]);
sigma=getopt(opts,'sigma',0.0);
B=getopt(opts,'massmat',speye(n));
deflMethod=0;
if ~isempty(seedSpace)
  deflMethod=2; % deflMethod 1: simple A-orthogonal deflation
                %            2: BNN (Tang et al paper)
                % so far nothing is working in terms of deflation.
                % This method is equivalent to balancing Neumann-Neumann,
                % and in the notation in Tang et. al. chooses as preconditioner
                % M1=P'M\P+Q, with E=(V'AV), Q=V(E\V'), P=I-AQ.
end

% diagonal and sub/super diagonal of the Lanczos method
lanD0=[];
lanD1=[];

%{
idrS=4;
idrsOpts=struct;
if ~isempty(seedSpace)
  idrsOpts.P=seedSpace;
  idrS=size(seedSpace,2);
end

[x,flag,relres,iter,resvec]=idrs(A-sigma*B,b,idrS,tol,maxIter,[],[],x0,idrsOpts);
return;
%}
%[x,flag,relres,iter,resvec]=gmres(A-sigma*B,b,50,tol,maxIter,@apply_op,[],x0,M);
%return;
%fprintf('CARP-CG tol: %4.2e\n',tol);

debug=getopt(opts,'debug',false);


if (debug)
  tol
  maxIter
end

nrm_ai2=nrms_ai2(A',sigma);
bnul=zeros(n,1);


V=[];
AV=[];
E=[];
if (deflMethod~=0)
  V=seedSpace;
  %AV= A*V - sigma*V;
  AV=V-dkswp(A,sigma,B,bnul,V,omega,nrm_ai2);
  E=V'*AV;
end  

nrm_b=norm(b);
nrm_r0=norm(A*x0-sigma*x0-b);
reltol2=tol*tol*nrm_b*nrm_b;

x=x0;
r=dkswp(A,sigma,B,b,x,omega,nrm_ai2)-x;

if (deflMethod==2)
  Vtr=V'*r;
  vr=E\Vtr;
  rtil=r-AV*(vr);
  z=apply_op(rtil,M);
  z=z-AV*(E\(V'*z));
  z=z+V*vr;
else
  z=apply_op(r,M);

end
p=z;

if (deflMethod==1)
  p = p - V*(AV'*p);
end

r2_new = r'*z;

alpha_old=0;
beta_old=0;

disp(sprintf('%d\t%e\t%e',0,sqrt(r2_new),sqrt(r2_new)/nrm_b));
for k=1:maxIter 
  if((mod(k,cfreq) == 0) | (correction_needed==true))   %self stabilizing carp-cg, cfreq = correction frequency    
  %if(0)                              
    q = p-dkswp(A,sigma,B,bnul,p,omega,nrm_ai2);
    r=dkswp(A,sigma,B,b,x,omega,nrm_ai2)-x;
    alpha = (r'*p)/(p'*q);
    x=x+alpha*p;

    if (mod(k-1,itcheck)==0)
      nrm_r = norm(A*x-sigma*x-b);
       if (mod(k-1,itprint)==0)
        disp(sprintf('%d\t%e\t%e\t%e\t%e',k,nrm_r,nrm_r/nrm_b,nrm_r-tol*nrm_r0,nrm_r0));
      end
      relres=nrm_r/nrm_b;
      resvec=[resvec,nrm_r];

      if (nrm_r<=tol*nrm_r0)
        save residuum_100p.txt resvec -ASCII;
        break;
      end
    end
    r=r-alpha*q;

    if (deflMethod==2)
      Vtr=V'*r;
      vr=E\Vtr;
      rtil=r-AV*(vr);
      z=apply_op(rtil,M);
      z=z-AV*(E\(V'*z));
      z=z+V*vr;
    else
      z=apply_op(r,M);
    end
    r2_old=r2_new;
    r2_new=r'*z;
    beta = -(r'*q)/(p'*q);
    p = r + beta*p;

    if (deflMethod==1)
      p = p - V*(AV'*z);
    end
    relres=sqrt(r2_old)/nrm_b;

   if (wantT)
      lanD0=[lanD0;1/alpha];
      lanD1=[lanD1;sqrt(beta)/alpha];
      if (k>1)
        lanD0(k) = lanD0(k)+beta_old/alpha_old;
      end
      alpha_old=alpha;
      beta_old=beta;
    end
  else                                                 %normal carp-cg                                                                 
    q=p-dkswp(A,sigma,B,bnul,p,omega,nrm_ai2);

    %provoke pseudo-random error
    if errormode
      q = elem_destr(q);     %endless loop
    end

    alpha = (r'*z)/(p'*q);
    x=x+alpha*p;

    if (mod(k-1,itcheck)==0)
      nrm_r = norm(A*x-sigma*x-b);
      if (mod(k-1,itprint)==0)
        disp(sprintf('%d\t%e\t%e\t%e',k,nrm_r,nrm_r/nrm_b,nrm_r-tol*nrm_r0));
      end
      relres=nrm_r/nrm_b;
      resvec=[resvec,nrm_r];

      if (nrm_r<tol*nrm_r0)
        save residuum_100p.txt resvec -ASCII;
        break;
      end
    end
    r=r-alpha*q;

    %provoke pseudo-random error
    
    if errormode
      for i = 1:10
        r = elem_destr(r);       %endless loop
      end
    end

    if (deflMethod==2)
      Vtr=V'*r;
      vr=E\Vtr;
      rtil=r-AV*(vr);
      z=apply_op(rtil,M);
      z=z-AV*(E\(V'*z));
      z=z+V*vr;
    else
      z=apply_op(r,M);
      %provoke pseudo-random error
      if errormode
      %for i = 1:100
        z = elem_destr(z);    %plenty more steps, slows the process a lot
      %end
      end
    end
    r2_old=r2_new;

    r2_new=r'*z;
    beta=r2_new/r2_old;
    p=z+beta*p;

    %provoke pseudo-random error
    if errormode
      p = elem_destr(p);    %plenty more time, D0 = 0, E0 = other value
    end

    if (deflMethod==1)
      p = p - V*(AV'*z);
    end
    relres=sqrt(r2_old)/nrm_b;
    %fprintf('\t%d\t%e\n',k,relres);
    %if (r2_old<reltol2) 
    %  break;
    %end
  if (wantT)
      lanD0=[lanD0;1/alpha];
      lanD1=[lanD1;sqrt(beta)/alpha];

      %provoke pseudo-random error
      if errormode
      %for i = 1:1000
         lanD1 = elem_destr(lanD1);    %plenty more time, E0 = 0
      %end
      end
      if (k>1)
        lanD0(k) = lanD0(k)+beta_old/alpha_old;

        %provoke pseudo-random error
        if errormode
        %for i = 1:1000
          lanD0 = elem_destr(lanD0);    %plenty more time, D0 = 0, E0 = other value
        %end
        end
      end
      alpha_old=alpha;
      beta_old=beta;
    end
  end
end
iter=k

if iter>=maxIter
  flag=1;
else
  flag=0;
end

if (wantT)
  Tlan = spdiags([lanD1,lanD0,[0;conj(lanD1(1:k-2))]],-1:1,k-1,k-1);
  disp(['D0=',num2str(lanD0(1))]);
  disp(['E0=',num2str(lanD1(1))]);
%Tlan=[lanD1, lanD0,[0;lanD1(1:k-1)]];
end
end
