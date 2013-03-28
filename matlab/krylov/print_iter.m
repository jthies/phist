function print_iter(iter,relres)

global printOpts;
debug=getopt(printOpts,'debug',false);
tab='\t';
fmt='';
if (isfield(printOpts,'indent'))
  for i=1:printOpts.indent
    fmt=[fmt,tab];
  end
end

fmt=[fmt,'%d'];

k=length(relres);
for i=1:k
  fmt=[fmt,tab,'%8.4e'];
end
if debug
  disp('=======================================');
end
disp(sprintf(fmt,iter,relres));
if debug
  disp('=======================================');
end

end
