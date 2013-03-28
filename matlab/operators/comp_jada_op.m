function M=comp_jada_op(A,sigma,B, V)
% computes X=(B-VV')(A-sigma*B)(B-VV')

n=size(A,1);
M.label='(I-vvT)(A-sB)(I-vvT)';

M.arg{1}=comp_proj_op(V);
M.arg{2}=comp_shifted_op(A,-sigma,B);
M.arg{3}=comp_proj_op(V);

M.compute=@comp_jada_op;
M.apply=@apply_sequence;


end
