function jd=CheckSortSchur(jd,sigma)

k=size(jd.Rs,1); 
if k==0 
  return
end

I=SortEig(diag(jd.Rs),sigma);

if ~min((1:k)'==I)
   [U,jd.Rs]=SwapSchur(eye(k),jd.Rs,I);
   jd.Qs=jd.Qs*U;
end

return
end
