function X=jada_op(Y, A,sigma,B, V)
% computes X=(B-VV')(A-sigma*B)(B-VV')

X = B*Y-(V'*Y)*V;
X = A*X - sigma*X;
X = B*X-(V'*X)*V;

end
