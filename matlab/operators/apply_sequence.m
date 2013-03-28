function X=apply_sequence(Y,varargin);

X=Y;
for i=1:nargin-1
  X=apply_op(X,varargin{i});
end

end
