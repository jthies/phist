function nrms=nrms_ai2(A,sigma)

global is_matlab is_octave

if length(is_matlab)==0
  is_matlab=~is_octave;
end

if (is_matlab)
  nrms=nrms_ai2_c(A,sigma);
else
  n=size(A,1);
  m=size(sigma,2);
  if (size(sigma,1)>1)
    error('nrms_ai2: sigma must be a row vector');
  end
  nrms=kron(-2*sigma,full(spdiags(A,0))+sigma*sigma)
  for i=1:n
    tmp=A(i,:)*A(i,:)';
    nrms(i,:)=nrms(i,:)+repmat(tmp,1,m);
  end
end



end
