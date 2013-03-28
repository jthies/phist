function M=comp_iluprec(A,M)

M.label='ilu(0)';
%opts.type='nofill';
opts.type='ilutp';
opts.droptol=0.001;
[M.arg{1},M.arg{2}]=ilu(A,opts);
M.apply=@apply_lu;
M.compute=@comp_iluprec;

end
