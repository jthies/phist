function print_eigs_iter(iter,m,res,val)

global printOpts;
tab='\t';
fmt='';
if (isfield(printOpts,'indent'))
  for i=1:printOpts.indent
    fmt=[fmt,tab];
  end
end

iscomplex=abs(imag(val))>sqrt(eps);

fmt=[fmt,'%d %d\t%12.5e '];

if iscomplex
  fmt=[fmt,' %+12.5e i'];
else
  fmt=[fmt,'\t'];
end

fmt=[fmt,'\t||r||=%16.8e'];
%disp('=================================================');
if iscomplex
  disp(sprintf(fmt,iter,m,real(val),imag(val),res));
else
  disp(sprintf(fmt,iter,m,real(val),res));
end
%disp('=================================================');

end
