function M=comp_AAt_op(A,sigma)
%function M=comp_AAt_op(A,sigma,omega)

n=size(A,1);
M.label='(A-sigma I)*(A-sigma I)^T';

M.arg{1}=A;
M.arg{2}=sigma;

M.compute=@comp_AAt_op;
M.apply=@apply_AAt;


end
