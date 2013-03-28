function [alpha,beta]=house(y)
%                                       
%function [alpha]=house(y)              
%                                       
% returns scalars alpha, beta such that 
% (I+beta*vv')y=alpha*e1, with          
% v=y - alpha*e1. y should be a vector. 
%                                       
% cf. Gutknecht, Schmelzer (2005),      
% alg. 5.                               
%                                       
ny=norm(y,2);
beta=0;
alpha=y(1);
if (ny~=0)
  if (y(1)~=0)
    f=y(1)/abs(y(1));
  else
    f=1;
  end
  alpha=-ny*f;
  beta=-1./(ny*(ny+abs(y(1))));
end

end
