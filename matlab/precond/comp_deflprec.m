function M=comp_deflprec(M0,T,U)
%function M=comp_deflprec(M0,T,U)
% deflation by skey projection, M*(I-U*(T\U'))
M.label=sprintf('deflation, %d vectors',size(U,2));
M.arg{1}=M0;
M.arg{2}=T;
M.arg{3}=U;
M.apply=@apply_deflprec;
end
