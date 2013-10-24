function A=make_testmat(nx,sym)

if (nargin<2)
  sym=false;
end

n=nx^3;

e=ones(nx,1);
A1=spdiags([-e 2*e -e],-1:1,nx,nx);
I1=speye(nx,nx);

A2=kron(I1,A1)+kron(A1,I1);
I2=speye(nx*nx,nx*nx);
A=kron(I1,A2)+kron(A1,I2);

if (sym)
  R=sprandsym(A);
else
  R=sprandn(A);
end
A=A+R;

end
