function M=comp_iterprec(A,M)

M.label='GMRES';
M.arg{1}=A;
if isfield(M,'opts')
  M.arg{2}=M.opts;
end
if isfield(M,'M')
  M.arg{3}=M;
else
  M.arg{3}=comp_idprec(A);
end

M.apply=@apply_iterprec;
M.compute=@comp_iterprec;
end
