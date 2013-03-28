function [Hcol,transj,rsjjp1]=updateQR(Hj,Hcol,trans,rsjjp1)
%                                                               
%function [Hcol,transj,rsjjp1]=updateQR(Hj,Hcol,trans,rsjjp1)   
%                                                               
% input: previous j x j-1 (j*k x (j-1)*k) tansformed            
%       H-matrix Hj,                                            
%        new column Hcol (dim. (j+1)*k x k)                     
%        transformation coefficients (Givens or HH) for         
%        previous j-1 rows of H                                 
%        right-hand side to be transformed along with H [j,j+1] 
%                                                        
% output Hcol and rs overwrite inputs, trans(:,j)=transj.
%                                                        
global printOpts;
debug=getopt(printOpts,'debug',false);
k=size(rsjjp1,2);
j=size(Hj,1)/k;
if (debug)
j
Hj
Hcol
rsjjp1
disp('Giv');
j
disp('length(trans)');
size(trans)
end

idx=@(jj) ((jj-1)*k+1):(jj*k);
jdx=idx(j);

% apply previous (j-1) transformations to columns j
if (k==1)
  for jj=1:size(trans,2) % apply Givens rotation
  %disp(sprintf('%d %f %f',jj,trans(1,jj),Hcol(jj)));
  %disp(sprintf('\t %f %f',trans(2,jj),Hcol(jj+1)));
    htmp = trans(1,jj)*Hcol(jj) + ...
           trans(2,jj)*Hcol(jj+1);
    Hcol(jj+1) = -conj(trans(2,jj))*Hcol(jj) + trans(1,jj)*Hcol(jj+1);
    Hcol(jj)   = htmp;
  end
else
  % Householder
  for jj=1:size(trans,2)
    jjdx=idx(jj);
    for kk=1:k
      rows=jjdx(kk):jjdx(end)+kk;
      beta=trans(1,jj,kk);
      v=Hj(rows,rows(1));
      v(1)=trans(2,jj,kk);
      Hcol(rows,:)=Hcol(rows,:) +beta*v*(v'*Hcol(rows,:));
    end
  end
end

% update QR factorization of H
if (k==1)
  % new Givens rotation for eliminating H(j+1,j)
  [cs,sn,htmp]=give(Hcol(j),Hcol(j+1));
  transj(1,1)=cs;
  transj(2,1)=sn;
  % eliminate H(j+1,j)
  Hcol(j) = htmp;
  Hcol(j+1) = 0.0;
  % apply to RHS
  tmp=cs*rsjjp1(idx(1));
  rsjjp1(idx(2)) = -conj(sn).*rsjjp1(idx(1));
  rsjjp1(idx(1))=tmp;
else
  % Householder transform to clear out the subdiagonal entries
  for kk=1:k
    rows=jdx(kk):jdx(end)+kk;
    cols=kk:k;
    v=Hcol(rows,kk);
    [alpha,beta]=house(v);
    transj(1,1,kk)=beta;
    v(1)=v(1)-alpha;
    transj(2,1,kk)=v(1);
    Hcol(rows(1),cols(1)) = alpha;
    Hcol(rows,cols(2:end)) = Hcol(rows,cols(2:end)) +beta*v*(v'*Hcol(rows,cols(2:end)));
    % apply to RHS
    rsjjp1(kk+[0:k],:) = rsjjp1(kk+[0:k],:) + ...
                beta*v*(v'*rsjjp1(kk+[0:k],:));
  end
end
if (debug)
disp('after transform');
Hcol
rsjjp1
disp('new Givens coeffs:');
transj.'
end
end
