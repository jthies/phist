function x=dkswp(A,b,x,omega,nrm_ai2)

% forward Kaczmarz sweep (CARP is Kaczmarz in parallel)
n=size(A,1);
for i=1:n
  scal=(b(i) - A(i,:)*x)/nrm_ai2(i);
  x = x + omega*scal*A(i,:)';
end

% ... and back ...
for i=n:-1:1
  scal=(b(i) - A(i,:)*x)/nrm_ai2(i);
  x = x + omega*scal*A(i,:)';
end

end
