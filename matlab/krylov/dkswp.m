function x=dkswp(A,sigma,B,b,x,omega,nrm_ai2)
%function x=dkswp(A,sigma,B,b,x,omega,nrm_ai2)
% Kaczmarz forward-backward sweep for system (A-sigma*I)x=b,

global nm_operations is_octave is_matlab

if (length(is_matlab)==0&&length(is_octave)==0)
  error('before calling dkswp, you must set is_matlab or is_octave');
else if (length(is_matlab)==0)
  is_matlab=~is_octave;
end
if isempty(nm_operations)
  nm_operations=0;
end

if (is_matlab)

% call mex function
x=dkswp2(A', sigma, B, b, x, omega, nrm_ai2);
return;

else

% slow reference implementation

n=size(A,1);
A=A-sigma*speye(n);


%tic;

% forward Kaczmarz sweep (CARP is Kaczmarz in parallel)
for i=1:n
  idx=find(A(i,:));
  scal=(A(i,idx)*x(idx)-b(i))/nrm_ai2(i);
  x(idx) = x(idx) - omega*scal*A(i,idx)';
end

% ... and back ...
for i=n:-1:1
  idx=find(A(i,:));
  scal=(A(i,idx)*x(idx)-b(i))/nrm_ai2(i);
  x(idx) = x(idx) - omega*scal*A(i,idx)';
end

%toc
%

end

nm_operations=nm_operations+1;

end
