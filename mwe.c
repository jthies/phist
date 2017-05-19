#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#include <omp.h>

/* minimal working example */

#include "esmac_types.h"
#include "esmac_generator.h"
#include "esmac_collection.h"

// split integer range [a...b-1] in n nearly equally sized pieces [ia...ib-1], for i=0,...,n-1
void split_range(esmac_idx_t a, esmac_idx_t b, esmac_idx_t n, esmac_idx_t i, esmac_idx_t *ia, esmac_idx_t *ib) {
  esmac_idx_t m = (b-a-1)/n + 1;
  esmac_idx_t d = n-(n*m -(b-a));
  if (i < d) {
    *ia = m*i + a;
    *ib = m*(i+1) + a;
  } else {
    *ia = m*d + (i-d)*(m-1) + a;
    *ib = m*d + (i-d+1)*(m-1) + a;
  }
}

int main(int argc, char *argv[]) {

  char *examplestr = argv[1];

  int info;
  esmac_generator_t * my_gen = esmac_collection_parse(examplestr, &info);

	printf("Example\n-------\n%s\n",my_gen->name);
	
	char * data_s = esmac_collection_example_data(my_gen);
      
	printf("\nParameters\n----------\n%s\n\n",data_s);
  
  esmac_idx_t nrow=0, n_nz=0;
  double elapsedt=0.0;

  #pragma omp parallel
  {
    double t1 = omp_get_wtime();
    
    esmac_generator_work_t * my_ws = esmac_generator_alloc(my_gen, &info);
 
    esmac_idx_t my_nrow = esmac_generator_query(my_gen,my_ws,"nrow");
    esmac_idx_t mrow = esmac_generator_query(my_gen,my_ws,"maxnzrow");

    esmac_idx_t *cind = malloc(mrow * sizeof *cind);
    double *val = malloc(mrow * sizeof *val);
  
    esmac_idx_t idx;
    esmac_idx_t ia,ib;
    // current thread will generate rows ia ... ib-1
    split_range(0,my_nrow, omp_get_num_threads(), omp_get_thread_num(), &ia, &ib);

    esmac_idx_t my_nz = 0;


    for (idx=ia;idx<ib;idx++) {
      int k = esmac_generator_row(my_gen, my_ws, idx, cind, val);
      my_nz += k;
    }

    free(cind);
    free(val);
    esmac_generator_free(my_ws);
        
    double t2 = omp_get_wtime();
    
    #pragma omp critical 
    {
      printf("========================================\nThread number: %d [of %d threads]\n generated %"ESMACPRIDX" rows with %"ESMACPRIDX" non-zeros\n in %12.3f seconds.\n",
            omp_get_thread_num(),omp_get_num_threads(), ib-ia, my_nz, t2-t1); 
    
      nrow += ib-ia;
      n_nz += my_nz;
      elapsedt = fmax(elapsedt,t2-t1);
    }
    
  }
    
  printf("\n========================================\nGenerated a total of %"ESMACPRIDX" rows with %"ESMACPRIDX" non-zeros\n in %12.3f seconds [%d rows/second]\n",
    nrow, n_nz, elapsedt, (int) ((double) nrow/elapsedt));  
	

  esmac_collection_example_free(my_gen);

  return 0;
}
