function [c, s, r] = give(f, g)
%function [cs, sn, r] = give( a, b )
% returns Givens rotation such that
% | c s ||f| |r|
% | _   || |=| |
% |-s c ||g| |0|
% This file uses the convention for the BLAS,
% cf. for instance Bindel et al. 2001). The  
% BLAS routine is called XROTG.

d=sqrt(abs(f).^2 + abs(g).^2);
c=zeros(size(f));
s=c; r=c;

for k=1:length(f)
   if ( g(k) == 0.0 ) % includes case f=0 and g=0
     cs(k) = 1.0;
     sn(k) = 0.0;
     r(k)  = f(k);
   elseif (f(k) == 0)
     c(k)=0;
     s(k)=sign(conj(g(k)));
     r(k)=abs(g(k));
   else
     c(k)=abs(f(k))/d(k);
     s(k)=sign(f(k))*conj(g(k))/d(k);
     r(k)=sign(f(k))*d(k);
   end
end

end
