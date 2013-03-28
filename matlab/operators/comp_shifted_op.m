function M=comp_shifted_op(A, sigma,B)
%function M=comp_shifted_op(A, sigma,B)
% constructs A+sigma*B
M.label='shifted op';
M.arg{1}=A;
M.arg{2}=sigma;
M.arg{3}=B;
M.apply=@apply_shifted_op;
M.compute=@comp_shifted_op;
end
