#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "esmac_cmk.h"


/* compare first index of pair of ints */
static int cmp_fst(  const void *a, const void *b ){ 
  int x,y;
  x=*(int *) a;
  y=*(int *) b;
  // prone to overflow - but not here.
  return (a-b);
}


int * esmac_cmk(const esmac_sparsemat_t *sm) {

  // matrix dimension
  int dim = sm->nr;

  // --- does not account for duplicates --- doesn't matter, either
  int maxrowlength= esmac_sparsemat_maxrowlength(sm);
  
  int * perm = malloc(dim * sizeof *perm);
  
  // work arrays
  int * vdeg   =malloc(dim * sizeof * vdeg);
  int * nodelist=malloc(dim * sizeof * nodelist);
  int * sortarray=malloc(2 * maxrowlength * sizeof * sortarray);
     
  int i, idx;
  
  // vertex degree
  for (i=0;i<dim;i++) {
    vdeg[i]=sm->rptr[i+1]-sm->rptr[i];
  }
      
    /*    
    if (!vdeg || !alist || !newlist || !newvd || !isdone || !sortarray) {
        printf("abort -> malloc failed!\n");
        exit(0);
    }
    */
   

  //create permutation
  
  //add  node with smallest vertex degree first
    
  int dummy=vdeg[0];
  idx=0;
  for (i=1;i<dim;i++) {
    if (vdeg[i]<dummy) {
			dummy=vdeg[i];
			idx=i;
		}
	}
  // printf("Node with smallest vertex degree: %d\n",idx);
    
  // perm not yet set
  for (i=0;i<dim;i++) {perm[i]=-1;}
    
  nodelist[0]=idx;
  perm[idx]=0;
  int ndone=1;
  
  for (idx=0;idx<dim;idx++) {
    int nidx=nodelist[idx];  // !from 0
    //construct list of *new* nodes, reachable from nidx
    int nnew=0;
    for (i=0;i<sm->rptr[nidx+1]-sm->rptr[nidx];i++) {
      if (perm[sm->cind[i+sm->rptr[nidx]]] < 0) {
        nodelist[ndone+nnew]=sm->cind[i+sm->rptr[nidx]];
        nnew++;
      }
    }
    if (nnew + ndone > dim) {
      printf("%s: problem with nnew + idx > dim\n",__func__);
      exit(EXIT_FAILURE);
    }   
    //sort new nodes in order of increasing vertex degree
  
    // we cheat a little.
    for (i=0;i<nnew;i++) {
		  sortarray[2*i  ]=vdeg[nodelist[ndone+i]];
		  sortarray[2*i+1]=nodelist[ndone+i];
    }
    qsort(sortarray,nnew, 2* sizeof * sortarray,cmp_fst);
    //and add them
    for (i=0;i<nnew;i++) {
      nodelist[ndone+i]=sortarray[2*i+1];
      perm[sortarray[2*i+1]]=ndone+i;
    }
       
    ndone=ndone+nnew;
  }
    
  free(vdeg);
	free(nodelist);
	free(sortarray);

  return perm;
    
}


esmac_sparsemat_t * esmac_sparsemat_permute(const esmac_sparsemat_t *sm, const int *perm) {

  if (sm->nr != sm->nc) {
    printf("%s: Require square matrix\n",__func__);
    exit(EXIT_FAILURE);
  }
  
  int dim = sm->nr;

  // inverse permutation
  int * invperm = malloc(sm->nr * sizeof *invperm);
  int i;
  for (i=0;i<dim;i++) {
    invperm[perm[i]]=i;
  }

/*
  int * dummy;
  dummy = perm;
  perm = invperm;
  invperm=perm;
*/

  esmac_sparsemat_t * sm_perm = esmac_sparsemat_alloc(sm->nr,sm->nc,sm->ne);

  esmac_idx_t n = 0;
  sm_perm->rptr[0]=0;
  esmac_idx_t idx;
  for (idx=0;idx<sm->nr;idx++) {
    esmac_idx_t orig_idx = invperm[idx];
    int k = sm->rptr[orig_idx+1] - sm->rptr[orig_idx];
    //beware about int <-> esmac_idx_t
    int i;
    for (i=0;i<k;i++) {
      sm_perm->cind[n+i]=perm[sm->cind[sm->rptr[orig_idx]+i]];
    }
    memcpy(&(sm_perm->val[n]), &(sm->val[sm->rptr[orig_idx]]) , k * sizeof *(sm->val));
    n=n+k;
    sm_perm->rptr[idx+1]=n;
  }
  if (n != sm->ne) {
    printf("%s: Error.",__func__);
    exit(EXIT_FAILURE);
  }
  sm_perm->ne=sm->ne;
  
  free(invperm);
  return sm_perm;
}

