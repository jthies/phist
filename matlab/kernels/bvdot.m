function h=bvdot(V,W)
%function h=bvdot(V,W)
% computes delta_ij*v_i'w_j
m=size(V,2);
h=zeros(1,m);

% parallel for
for i=1:m
  h(i) = V(:,i)'*W(:,i);
end

% in an MPI context, the communication for the dot-products can be bundled.
% MPI_Allreduce(h,+)

end
