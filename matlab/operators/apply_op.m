function x=apply_op(y,op)
%function x=apply_op(y,op)
% applies some generally defined operator of the form
% x=op.fun(y,op.arg{:}). If op is not a struct, we try
% to return x=op*y instead.

global count_MVM;

if isempty(y)
  x=y;
  return;
end

if (isstruct(op))
%  disp(['apply ',op.label]);
  x=op.apply(y,op.arg{:});
else
  x=op*y;  
  if issparse(op)
    % check if this is the identity matrix
    % (this is of course terribly expensive
    % and just for getting the counts right,
    % in a practical implementation we'ld do
    % it smarter)
    is_I=false;
    if (norm(op-speye(size(op)),inf)==0)
      is_I=true;
    end
    if (~is_I)
      if isempty(count_MVM)
        count_MVM=0;
      end
      count_MVM=count_MVM+1;
    end
  end
end

end
