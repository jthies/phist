#include <stdlib.h>

#include "esmac_lut.h"

    

esmac_idx_t esmac_lut_construct(int ineq, int n, int m, int s, esmac_idx_t **lut) {
	int i,j,k;
	
  *lut = malloc(n*(s+1) * sizeof (esmac_idx_t));
  
  esmac_idx_t *cnt = *lut;
	
	for (i=0;i<n*(s+1);i++)
	{ 
		cnt[i]=0;
	}
	
	for (i=0;i<=m && i<=s;i++)
	{
		cnt[i]=1;
	}
	for (j=2;j<=n;j++)
	{
		for (k=0;k<=s;k++)
		{
			for (i=0;i<=m && i<=k;i++)
			{
				cnt[(s+1)*(j-1)+k]=cnt[(s+1)*(j-1)+k]+cnt[(s+1)*(j-2)+(k-i)];
			}
		}
	}
	/* total number of states */
	esmac_idx_t ntot;
	if (ineq)
	{
		ntot=0;
		for (i=0;i<=s;i++)
		{
			ntot=ntot+cnt[(n-1)*(s+1)+i];
		}
	}
	else
	{
		ntot=cnt[n*(s+1)-1];
	}
	
	return ntot;
}


// TODO: include check on range of x[i] \in 0,..,M-1
// this routine (c/s)hould be optimized further
esmac_idx_t esmac_lut_encode(int ineq, int n, int s, esmac_idx_t *cnt, const int *x) {
	int i,j;
	int sum;
	esmac_idx_t idx;
	
	idx=0;
	sum=x[n-1];
	
	for (i=2;i<=n;i++)
	{
		sum=sum+x[n-i];
		for (j=0;j<x[n-i];j++)
		{
			idx=idx+cnt[(s+1)*(i-2)+(sum-j)];
		}
	}
	
	if (ineq)
	{
		for (j=0;j<sum;j++)
		{
			idx=idx+cnt[(n-1)*(s+1)+j];
		}
	}

	return idx;
}


int esmac_lut_decode(int ineq, int n, int s, esmac_idx_t *cnt, esmac_idx_t idx, int *x) {
	int i,nn;
	int sum;
	
	if (ineq)
	{
		sum=0;
		while (idx>=cnt[(n-1)*(s+1)+sum])
		{
			idx=idx-cnt[(n-1)*(s+1)+sum];
			sum++;
		}
	}
	else
	{
		sum=s;
	}
		
	for (nn=1;nn<=n;nn++)
	{
		if (nn==n)
		{
			x[nn-1]=sum;
		}
		else
		{
			x[nn-1]=0;
			i=sum;
			while (idx >= cnt[(s+1)*(n-nn-1)+i])
			{
				x[nn-1]++;
				idx=idx-cnt[(s+1)*(n-nn-1)+i];
				i=i-1;
			}
			sum=sum-x[nn-1];
		}
	}
	
	return 0;
}


esmac_idx_t esmac_lut_onezero_construct(int n, int s, esmac_idx_t **lut) {
  /*
  int *lut = malloc(n*(s+1) * sizeof (esmac_idx_t));
  int *my_lut;
  esmac_lut_construct(0, n, 1, s, &my_lut);
  lut[0]=n;
  lut[1]=s;
  int i;
  for (i=0;i<  ;i++) {
    lut[i+2] = my_lut[(s+1)* i + ...  ];
  free(my_lut);
  */
  return esmac_lut_construct(0, n, 1, s, lut);
}

esmac_idx_t esmac_lut_onezero_encode(int n, int s, esmac_idx_t *cnt, const int *x) {
  int i;
	int sum;
	esmac_idx_t idx;
	
	idx=0;
	sum=x[n-1];
	
	for (i=2;i<=n;i++) {
    if (x[n-i]) {
      sum++;
			idx=idx+cnt[(s+1)*(i-2)+sum];
      if (sum==s) {break;}
		}
	}

	return idx;
  
}
	  
int esmac_lut_onezero_decode(int n, int s, esmac_idx_t *cnt, esmac_idx_t idx, int *x) {
  int nn;
	int sum;
	
  sum=s;
		
	for (nn=0;nn<n-1;nn++) {
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
	
	return 0;
  
}

