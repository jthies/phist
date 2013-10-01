function jd=CheckForNullSpace(jd,nselect,tau,tol,W,E,Rv);
% V,W orthonormal, A*V-tau*V=W*R+jd.Qs'*E

  k=size(jd.Rs,1); 
  j=size(V,2);

  [jd.W,jd.R]=qr(jd.W,0); 
  E=E/Rv; 
  R=R/Rv; 
  jd.M=jd.W'*jd.V;
  %%% not accurate enough M=Rw'\(M/Rv);

  if k>=nselect
    return
  end

  CHECK=1; 
  l=k;

  [S,T,Z,Q]=qz(R,M); 
  Z=Z';
  while CHECK
    I=SortEigPairVar(S,T,2); 
    [Q,Z,S,T]=SwapQZ(Q,Z,S,T,I(1));
    s=abs(S(1,1)); 
    t=min(abs(T(1,1)),1); 
    CHECK=(s*sqrt(1-t*t)<tol);
    if CHECK
      jd.V=jd.V*Q; 
      jd.W=jd.W*Z; 
      jd.E=jd.E*Q;
      u=jd.V(:,1);
      [r,a]=RepGS(u,(jd.W(:,1)-T(1,1)'*u)*S(1,1),0);
      jd.Qs=[jd.Qs,u]; 
      t=(T(1,1)'-a/S(1,1))*S(1,:);
      Rs=[Rs,jd.E(:,1); zeros(1,k),tau+t(1,1)]; 
      k=k+1;
      J=[2:j]; j=j-1;
      jd.V=jd.V(:,J); 
      jd.W=jd.W(:,J); 
      jd.E=[jd.E(:,J);
      t(1,J)];  
      s=S(1,J)/S(1,1);
      jd.R=S(J,J); 
      jd.M=T(J,J); 
      Q=eye(j); 
      Z=eye(j);

      s=s/R; 
      nrs=norm(r)*norm(s);
      if nrs>=tol,
        jd.W=jd.W+r*s; 
        jd.M=jd.M+s'*(r'*V);
        if (nrs^2>eps) 
          [jd.W,R0]=qr(W,0); 
          jd.R=R0*jd.R; 
          jd.M=R0'\jd.M; 
        end
      end

      S=jd.R; T=jd.M;

      CHECK=(k<nselect & j>0);
    end
  end

end
