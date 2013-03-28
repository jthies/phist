function Y=apply_iterprec(X,A,opts,M)

global printOpts;
if ~isfield(printOpts,'indent')
  printOpts.indent=0;
end

printOpts.indent=printOpts.indent+1;
Y=bgmres(A,X,[],opts,M);
printOpts.indent=printOpts.indent-1;

end
