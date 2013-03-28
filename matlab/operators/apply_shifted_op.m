function X=apply_shifted_op(Y, A, sigma, B)
% computes X=(A+sigma*B)Y
X = apply_op(Y,A);
  
if (sigma)
  X=X+sigma*apply_op(Y,B);
end

end
