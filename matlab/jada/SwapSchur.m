function [Q,S]=SwapSchur(Q,S,I)
% [Q,S]=SwapSchur(QQ,SS,P)
%    QQ and SS are square matrices of size K by K
%    P is the first part of a permutation of (1:K)'.
%
%    If    M = QQ*SS*QQ'  and  QQ'*QQ = EYE(K), SS upper triangular
%    then  M*Q = Q*S      with   Q'*Q = EYE(K),  S upper triangular
%    and   D(1:LENGTH(P))=DD(P) where D=diag(S), DD=diag(SS)
%
%    Computations uses Givens rotations.

  kk=min(length(I),size(S,1)-1);
  j=1; while (j<=kk & j==I(j)), j=j+1; end;
  while j<=kk
    i=I(j);
    for k=i-1:-1:j
      q = [S(k,k)-S(k+1,k+1),S(k,k+1)];
      if q(1) ~= 0
        q = q/norm(q);
        G = [[q(2);-q(1)],q'];
        J = [k,k+1];
        Q(:,J) = Q(:,J)*G;
        S(:,J) = S(:,J)*G;
        S(J,:) = G'*S(J,:);
      end
      S(k+1,k) = 0;
    end
    I=I+(I<i);
    j=j+1; while (j<=kk & j==I(j)), j=j+1; end
  end

return
end
