function value=getopt(opts,name,default)
% function value=getopt(opts,name,default)
value=default;
if isfield(opts,name)
  value=getfield(opts,name);
end
end
