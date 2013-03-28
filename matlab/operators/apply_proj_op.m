function X=apply_proj_op(Y,V)

% computes X=(I-VV')Y
X = Y-V*(V'*Y);

end
