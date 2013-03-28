function M=comp_proj_op(V)

M.label='projection op';
M.arg{1} = V;
M.apply=@apply_proj_op;
M.compute=@comp_proj_op;

end
