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
    if isempty(count_MVM)
      count_MVM=0;
    end
    count_MVM=count_MVM+1;
  end
end

end
