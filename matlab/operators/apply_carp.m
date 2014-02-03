function X=apply_carp(Y,A,sigma,omega,diagA,nrm_ai2)

X=zeros(size(Y));

% forward Kaczmarz sweep (CARP is Kaczmarz in serial)
n=size(A,1);
for i=1:n
  % row indices
  scal=(Y(i) - A(i,:)*X)/nrm_ai2(i);
  X = X + omega*scal*A(i,:)';
end

% ... and back ...
for i=n:-1:1
  % row indices
  scal=(Y(i) - A(i,:)*X)/nrm_ai2(i);
  X = X + omega*scal*A(i,:)';
end

% I - (the whole thing)
X=Y-X;

end
