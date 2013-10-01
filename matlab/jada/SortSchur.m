%=======================================================================
%============== Sorts Schur form =======================================
%=======================================================================
function [Q,S]=SortSchur(A,sigma,gamma,kk)
%[Q,S]=SortSchur(A,sigma)
%  A*Q=Q*S with diag(S) in order prescribed by sigma.
%  If sigma is a scalar then with increasing distance from sigma.
%  If sigma is string then according to string
%  ('LM' with decreasing modulus, etc)
%
%[Q,S]=SortSchur(A,sigma,gamma,kk)
%  if gamma==0, sorts only for the leading element
%  else, sorts for the kk leading elements

  l=size(A,1);
  if l<2, Q=1;S=A; return,
  elseif nargin==2, kk=l-1;
  elseif gamma, kk=min(kk,l-1);
  else, kk=1; sigma=sigma(1,:); end

%%%------ compute schur form -------------
  [Q,S]=schur(A); %% A*Q=Q*S, Q'*Q=eye(size(A));
%%% transform real schur form to complex schur form
  if norm(tril(S,-1),1)>0, [Q,S]=rsf2csf(Q,S); end

%%%------ find order eigenvalues ---------------
  I = SortEig(diag(S),sigma);

%%%------ reorder schur form ----------------
  [Q,S] = SwapSchur(Q,S,I(1:kk));

return
end
