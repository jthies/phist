#include <stdlib.h>

#include "scamac_safeint.h"
#include "scamac_lut.h"



ScamacErrorCode scamac_lut_construct(bool ineq, int n, int m, int s, ScamacIdx * ns, ScamacIdx ** lut) {
  if (!ns) {
    return SCAMAC_ENULL;
  }
  if ( (n<=0) || (m<0) || (s<0) ) {
    return SCAMAC_ERANGE;
  }
  if ( (n>SCAMACHUGEINT) || (m>SCAMACHUGEINT) || (s>SCAMACHUGEINT) ) {
    return SCAMAC_EHUGEINT;
  }
  if ( (m==0) || (s==0) ) {
    *ns = 1;
    if (lut) {
      *lut = NULL;
    }
    return SCAMAC_EOK;
  }


  ScamacIdx nspo; // = n*(s+1)
  nspo = scamac_safe_mult(n,scamac_safe_add(s,1));
  if (nspo<0) {
    return SCAMAC_EOVERFLOW;
  }
  if (nspo> SCAMACHUGEINT ) {
    return SCAMAC_EHUGECOMP;  // avoid excessively large LUT
  }



  ScamacIdx *cnt = malloc( n*(s+1) * sizeof *cnt );
  if (cnt == NULL) {
    return SCAMAC_EMALLOCFAIL;
  }

  int i,j,k;

  for (i=0; i<n*(s+1); i++) {
    cnt[i]=0;
  }

  for (i=0; i<=m && i<=s; i++) {
    cnt[i]=1;
  }

  for (j=2; j<=n; j++) {
    for (k=0; k<=s; k++) {
      for (i=0; i<=m && i<=k; i++) {
        // cnt[(s+1)*(j-1)+k]=cnt[(s+1)*(j-1)+k]+cnt[(s+1)*(j-2)+(k-i)];
        cnt[(s+1)*(j-1)+k]=scamac_safe_add(cnt[(s+1)*(j-1)+k],cnt[(s+1)*(j-2)+(k-i)]);
        if (cnt[(s+1)*(j-1)+k] < 0) {
          free(cnt);
          return SCAMAC_EOVERFLOW;
        }
      }
    }
  }

  ScamacIdx ntot; // total number of states

  if (ineq) {
    ntot=0;
    for (i=0; i<=s; i++) {
      // ntot=ntot+cnt[(n-1)*(s+1)+i];
      ntot=scamac_safe_add(ntot,cnt[(n-1)*(s+1)+i]);
      if (ntot < 0) {
        free(cnt);
        return SCAMAC_EOVERFLOW;
      }
    }
  } else {
    ntot=cnt[n*(s+1)-1];
  }

  *ns = ntot;
  if (lut) {
    *lut=cnt;
  } else {
    free(cnt);
  }

  return SCAMAC_EOK;
}


// TODO: include check on range of x[i] \in 0,..,M-1
// this routine (c/s)hould be optimized further
ScamacIdx scamac_lut_encode(bool ineq, int n, int s, ScamacIdx * cnt, const int * x) {
  int i,j;
  int sum;
  ScamacIdx idx;

  idx=0;
  sum=x[n-1];

  for (i=2; i<=n; i++) {
    sum=sum+x[n-i];
    for (j=0; j<x[n-i]; j++) {
      idx=idx+cnt[(s+1)*(i-2)+(sum-j)];
    }
  }

  if (ineq) {
    for (j=0; j<sum; j++) {
      idx=idx+cnt[(n-1)*(s+1)+j];
    }
  }

  return idx;
}


void scamac_lut_decode(bool ineq, int n, int s, ScamacIdx * cnt, ScamacIdx idx, int * x) {
  int i,nn;
  int sum;

  if (ineq) {
    sum=0;
    while (idx>=cnt[(n-1)*(s+1)+sum]) {
      idx=idx-cnt[(n-1)*(s+1)+sum];
      sum++;
    }
  } else {
    sum=s;
  }

  for (nn=1; nn<=n; nn++) {
    if (nn==n) {
      x[nn-1]=sum;
    } else {
      x[nn-1]=0;
      i=sum;
      while (idx >= cnt[(s+1)*(n-nn-1)+i]) {
        x[nn-1]++;
        idx=idx-cnt[(s+1)*(n-nn-1)+i];
        i=i-1;
      }
      sum=sum-x[nn-1];
    }
  }

}


ScamacErrorCode scamac_lut_onezero_construct(int n, int s, ScamacIdx * ns, ScamacIdx ** lut) {
  return scamac_lut_construct(false, n, 1, s, ns, lut);
}

ScamacIdx scamac_lut_onezero_encode(int n, int s, ScamacIdx * cnt, const int * x) {
  int i;
  int sum;
  ScamacIdx idx;

  idx=0;
  sum=x[n-1];

  for (i=2; i<=n; i++) {
    if (x[n-i]) {
      sum++;
      idx=idx+cnt[(s+1)*(i-2)+sum];
      if (sum==s) {
        break;
      }
    }
  }

  return idx;

}

void scamac_lut_onezero_decode(int n, int s, ScamacIdx * cnt, ScamacIdx idx, int * x) {
  int nn;
  int sum;

  sum=s;

  for (nn=0; nn<n-1; nn++) {
    if (idx >= cnt[(s+1)*(n-nn-2)+sum]) {
      x[nn]=1;
      idx=idx-cnt[(s+1)*(n-nn-2)+sum];
      sum--;
      /*
      if (!sum) {
        int i;
        for (i=nn+1;i<n;i++) { x[i]=0; }
        return 0;
      }
      */
    } else {
      x[nn]=0;
    }
  }

  x[n-1]=sum;

}

