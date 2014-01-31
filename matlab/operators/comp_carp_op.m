function M=comp_carp_op(A,sigma,omega)
%function M=comp_carp_op(A,sigma,omega)

n=size(A,1);
M.label=sprintf('I-DKSWP(A-sigmaI,%4.2f)',omega);

M.arg{1}=A;
M.arg{2}=sigma;
M.arg{3}=omega;
M.arg{4}=full(diag(A));

nrm_ai2=zeros(n,1);

for i=1:n
  ai=A(i,:);
  ai(1,i)=ai(1,i)-sigma;
  nrm_ai2(i) = ai*ai';
end

M.arg{5}=nrm_ai2;

M.compute=@comp_carp_op;
M.apply=@apply_carp;


end
