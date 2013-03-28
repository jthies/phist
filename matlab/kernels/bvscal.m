function W=bvscal(V,s)
% function W=bvscal(V,s)
% computes W(:,i) = V(:,i)*s(i) for all i=1:size(V,2).

W=V;

for i=1:size(V,2)
  W(:,i) = W(:,i)*s(i);
end  

end
