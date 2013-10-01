function [v,y]=RepGS(V,v,gamma)
% [v,y]=REP_GS(V,w)
% If V orthonormal then [V,v] orthonormal and w=[V,v]*y;
% If size(V,2)=size(V,1) then w=V*y;
%
% The orthonormalisation uses repeated Gram-Schmidt
% with the Daniel-Gragg-Kaufman-Stewart (DGKS) criterion.
%
% [v,y]=REP_GS(V,w,GAMMA)
% GAMMA=1 (default) same as [v,y]=REP_GS(V,w)
% GAMMA=0, V'*v=zeros(size(V,2)) and  w = V*y+v (v is not normalized).


% coded by Gerard Sleijpen, August 28, 1998

if nargin < 3, gamma=1; end

[n,d]=size(V);

if size(v,2)==0, y=zeros(d,0); return, end

nr_o=norm(v); nr=eps*nr_o; y=zeros(d,1);
if d==0
  if gamma, v=v/nr_o; y=nr_o; else, y=zeros(0,1); end, return
end

y=V'*v; v=v-V*y; nr_n=norm(v); ort=0;

while (nr_n<0.5*nr_o & nr_n > nr)
  s=V'*v; v=v-V*s; y=y+s;
  nr_o=nr_n; nr_n=norm(v);     ort=ort+1;
end

if nr_n <= nr, if ort>2, disp(' dependence! '), end
  if gamma  % and size allows, expand with a random vector
    if d<n, v=RepGS(V,rand(n,1)); y=[y;0]; else, v=zeros(n,0); end
  else, v=0*v; end
elseif gamma, v=v/nr_n; y=[y;nr_n]; end

end
