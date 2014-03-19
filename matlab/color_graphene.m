function ord=color_graphene(nx,ny)

n=nx*ny; % matrix size
nc=4; % number of colors needed
nperc=n/nc;

members1=zeros(nperc,1);
members2=zeros(nperc,1);
members3=zeros(nperc,1);
members4=zeros(nperc,1);

gid = @(i,j) j*nx+i+1;
%gid = @(i,j) i*ny+j+1;

% distance-2 coloring with 4 colors
i0=0;
pos=ones(1,nc);

for j=0:2:ny-1
  for i=i0:2:nx-1
    members1(pos(1))=gid(i,j);
    pos(1)=pos(1)+1;
  end
  i0=1-i0;
  for i=i0:2:nx-1
    members2(pos(2))=gid(i,j);
    pos(2)=pos(2)+1;
  end
end
i0=0;
for j=1:2:ny-1
  for i=i0:2:nx-1
    members3(pos(3))=gid(i,j);
    pos(3)=pos(3)+1;
  end
  i0=1-i0;
  for i=i0:2:nx-1
    members4(pos(4))=gid(i,j);
    pos(4)=pos(4)+1;
  end
end

ord=[members1;
     members2;
     members3;
     members4];
