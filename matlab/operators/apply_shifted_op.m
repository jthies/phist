function X=apply_shifted_op(Y, A, sigma, B)
% computes X=(A+sigma*B)Y
X = apply_op(Y,A);
  
if (sigma~=0)
  X=X+sigma*apply_op(Y,B);
end

end
