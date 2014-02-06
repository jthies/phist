function x=dkswp(A,b,x,omega,nrm_ai2)
%function x=dkswp(A,b,x,omega,nrm_ai2)
% Kaczmarz forward-backward sweep for system Ax=b,

global nm_operations

% forward Kaczmarz sweep (CARP is Kaczmarz in parallel)
n=size(A,1);

for i=1:n
  scal=(A(i,:)*x-b(i))/nrm_ai2(i);
  x = x - omega*scal*A(i,:)';
end

% ... and back ...
for i=n:-1:1
  scal=(A(i,:)*x-b(i))/nrm_ai2(i);
  x = x - omega*scal*A(i,:)';
end

if isempty(nm_operations)
  nm_operations=0;
end
nm_operations=nm_operations+1;

end
