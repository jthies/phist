function X=jada_op(Y, A,sigma,B, V)
% computes X=(B-VV')(A-sigma*B)

X = A*Y - sigma*B*Y;
X = B*X-(V'*X)*V;

end
