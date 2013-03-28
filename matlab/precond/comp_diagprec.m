function M=comp_diagprec(A,varargin)
%function M=comp_diagprec(A)
%function M=comp_diagprec(A,M)


n=size(A,1);
M.label='I';
M.arg{1}=spdiags(1./spdiag(A,0),0,n,n);
M.apply=@apply_op;
M.compute=@comp_diagprec;

end
