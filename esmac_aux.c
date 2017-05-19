#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_aux.h"

int * esmac_counter_alloc(int n) {
  return malloc(n * sizeof (int));
}


void esmac_counter_reset(int n, int *c) {
  int i;
  for (i=0;i<n;i++) {
    c[i]=0;
  }
}

// increase counter by 1, such that c[i]<=maxc[i]
// returns 0 if counter reached last state
int esmac_counter_step(int n, const int *maxc, int *c) {
  int i = 0;
  while (i<n) {
    if (c[i]<maxc[i]) {
      c[i]=c[i]+1;
      int j;
      for (j=0;j<i;j++) {
	c[j]=0;
      }
      return i+1;
    }
    i=i+1;
  }
  return 0;
}

int esmac_increase_n_somewhat(int n) {
  return ((n >> 4) + (n>>6) + 2) << 4;
}


static int compare_longint(const void *a, const void *b) {
  esmac_idx_t mya = *((const esmac_idx_t *) a);
  esmac_idx_t myb = *((const esmac_idx_t *) b);
  if (mya < myb) {
    return -1;
  } else if (mya > myb) {
    return 1;
  } else {
    return 0;
  }
}

int esmac_sort_purge_array(int n, esmac_idx_t *arr) {
  if (n<=1) { return n;}

  // sort
  qsort(arr,n,sizeof *arr,compare_longint);

  // eliminate doublets
  
  int i,j;    
  i=0; 
  for (j=1;j<n;j++) {
    if (arr[i] != arr[j]) {
      i++;
      if (i<j) {
	arr[i] = arr[j];
      }
    }
  }

  return i+1;
	  
} 
