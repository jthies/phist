function M=comp_idprec(A,M)
%function M=comp_diagprec(A)
%function M=comp_diagprec(A,M)

M.label='identity';

n=size(A,1);
M.arg{1}=speye(n,n);
M.apply=@apply_op;
M.compute=@comp_idprec;

end
