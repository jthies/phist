function I=SortEig(t,sigma);
%I=SortEig(T,SIGMA) sorts the indices of T.
%
% T is a vector of scalars,
% SIGMA is a string or a vector of scalars.
% I is a permutation of (1:LENGTH(T))' such that:
%   if SIGMA is a vector of scalars then
%   for K=1,2,...,LENGTH(T) with KK = MIN(K,SIZE(SIGMA,1))
%      ABS( T(I(K))-SIGMA(KK) ) <= ABS( T(I(J))-SIGMA(KK) )
%      SIGMA(kk)=INF: ABS( T(I(K)) ) >= ABS( T(I(J)) )
%         for all J >= K

if ischar(sigma)
  switch sigma
    case 'LM'
      [s,I]=sort(-abs(t));
    case 'SM'
      [s,I]=sort(abs(t));
    case 'LR';
      [s,I]=sort(-real(t));
    case 'SR';
      [s,I]=sort(real(t));
    case 'BE';
      [s,I]=sort(real(t)); I=twistdim(I,1);
  end
else

  [s,I]=sort(abs(t-sigma(1,1)));
  ll=min(size(sigma,1),size(t,1)-1);
  for j=2:ll
    if sigma(j,1)==inf
      [s,J]=sort(abs(t(I(j:end)))); J=flipdim(J,1);
    else
      [s,J]=sort(abs(t(I(j:end))-sigma(j,1)));
    end
    I=[I(1:j-1);I(J+j-1)];
  end

end

return;
end
