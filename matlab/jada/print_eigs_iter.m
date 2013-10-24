function print_eigs_iter(iter,m,res,val)

global printOpts;
global count_MVM;
tab='\t';
fmt='';
if (isfield(printOpts,'indent'))
  for i=1:printOpts.indent
    fmt=[fmt,tab];
  end
end

iscomplex=true;
%iscomplex=abs(imag(val))>sqrt(eps);

fmt=[fmt,'it=%d m=%d #MV=%d\tth=%12.5e '];

if iscomplex
  fmt=[fmt,' %+12.5e i'];
else
  fmt=[fmt,'\t'];
end

fmt=[fmt,'\t||r||=%16.8e'];
%disp('=================================================');
if iscomplex
  disp(sprintf(fmt,iter,m,count_MVM,real(val),imag(val),res));
else
  disp(sprintf(fmt,iter,m,count_MVM,real(val),res));
end
%disp('=================================================');

end
