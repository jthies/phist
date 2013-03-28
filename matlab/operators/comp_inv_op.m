function M=comp_inv_op(A)
%function M=comp_inv_op(A)
% inv_op simply gives y=A\x
M.label='inv op';
M.arg{1}=A;
M.apply=@apply_inv_op;
M.compute=@comp_inv_op;
end
