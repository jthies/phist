n=25;
key='no key';

% all zero matrix

name=['spzero',int2str(n)];
description='zero matrix';

A=sparse(n,n);

write_all_formats(name,description,A,A);

% this fails, need to write the matrix manually
%dm2hb([name,'.rua'],A,description,key,'RUA',8,2);

% identity matrix

name=['speye',int2str(n)];
description='identity matrix';

A=speye(n,n);

write_all_formats(name,description,A,A);

% random matrices with and without zeros on diagonal, and row-sum 1

name=['sprandn_nodiag',int2str(n)];
description='random matrix with row sum 1 and some zeros on diagonal';

A=sprandn(n,n,0.25);
A(4,4)=randn; 
C=sprandn(A)+sprandn(A)*i;

rs=sum(A,2);
D=spdiags(1./rs,0,n,n);
A = D*A;

rs=sum(C,2);
D=spdiags(1./rs,0,n,n);
C = D*C;


% random matrices with and with some zeros on diagonal, and row-sum 1

name=['sprandn',int2str(n)];
description='random matrix with row sum 1 and nonzero diagonal';

A=sprandn(n,n,0.3);
A=sprandn(A+speye(n));

C=sprandn(A)+sprandn(A)*i;

rs=sum(A,2);
D=spdiags(1./rs,0,n,n);
A = D*A;

rs=sum(C,2);
D=spdiags(1./rs,0,n,n);
C = D*C;


write_all_formats(name,description,A,C);
