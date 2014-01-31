function X=apply_AAt(Y,A,sigma)

% TROET - unpreconditioned normal equations operator AA'
Ytmp=A'*Y-sigma*Y;
X = A*Ytmp-sigma*Ytmp;

end
