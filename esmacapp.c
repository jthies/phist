#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "esmac_config.h"
#include "esmac_generator.h"
#include "esmac_sparsemat.h"
#include "esmac_statistics.h"
#include "esmac_collection.h"

#ifdef ESMAC_USE_PNG
#include "esmac_plot.h"
#endif
#ifdef ESMAC_USE_LAPACK
#include "esmac_lanczos.h"
#endif
#ifdef ESMAC_USE_LAPACK
#include "esmac_spectrum.h"
#endif

#include "esmac_sparsemat_io.h"



const char * APPNAME = "esmacapp";

void print_help() {
  printf("%s - Exploring the ESMAC (ESsex MAtrix Collection)\n\n",APPNAME);
  printf("SYNOPSIS\n  %s command(s) [options] [EXAMPLE] [PARAMETERS]\n\n",APPNAME);
  printf("Example usage:\n");
  printf("  %s list\n",APPNAME);
  printf("  %s info EXAMPLE\n",APPNAME);
  printf("  %s query EXAMPLE [PARAMETERS]\n",APPNAME);
  printf("  %s stat EXAMPLE [PARAMETERS]\n",APPNAME);
  #ifdef ESMAC_USE_PNG
  printf("  %s plot EXAMPLE [PARAMETERS]\n",APPNAME);
  #endif
  #ifdef ESMAC_USE_LAPACK
  printf("  %s lanczos EXAMPLE [PARAMETERS]\n",APPNAME);
  #endif
  #ifdef ESMAC_USE_LAPACK
  printf("  %s spectrum EXAMPLE [PARAMETERS]\n",APPNAME);
  #endif
  printf("  %s output [mm|hb|mat|ghm] EXAMPLE [PARAMETERS]\n",APPNAME);
  printf("\nOptions:\n");
  printf("  -o, --outfile PREFIX\n");
  printf("  -q, --quiet\n");
  printf("  -p, --progress\n");
  printf("  --force\n");
  printf("  --maxdim DIM\n");
}

int main(int argc, char *argv[]) {

  esmac_generator_t * my_gen;

  int info;
  
  if (argc<2) {

    print_help();

  } else {
    
    if (!strcmp(argv[1],"list")) {
      
      char * s = esmac_collection_list_examples();
      
      printf("List of examples:\n\n%s\n\n",s);
    } else if (!strcmp(argv[1],"info")) {
      if (argc<3) {
        printf("\"info\" requires example\n");
        exit(EXIT_FAILURE);
      }
      char * sp = esmac_collection_list_parameters(argv[2], &info);
      if (info) {
        printf("unknown example \"%s\"\n",argv[2]);
        exit(EXIT_FAILURE);
      }
      printf("=== Example: %s ===\n\n=== parameters ===\n",argv[2]);
      printf("%s\n",sp);
    } else if (!strcmp(argv[1],"query")) {
      if (argc<3) {
        printf("\"query\" requires example\n");
        exit(EXIT_FAILURE);
      }
      my_gen = esmac_collection_parse(argv[2], &info);
      if (info) {
        printf("unknown example or parameters \"%s\" (try \"list\" or \"info\")\n",argv[2]);
        exit(EXIT_FAILURE);
      }
      esmac_generator_work_t * my_work = esmac_generator_alloc(my_gen, &info);
      char * data_s = esmac_collection_example_data(my_gen);
      printf("=== Example: %s ===\n\n=== parameters ===\n",argv[2]);
      printf("%s\n",data_s);
      printf("dimension nrow x ncol\n%"ESMACPRIDX" %"ESMACPRIDX"\n",esmac_generator_query(my_gen,my_work,"nrow"),esmac_generator_query(my_gen,my_work,"ncol"));
      free(data_s);
      esmac_generator_free(my_work);
      esmac_collection_example_free(my_gen);
    } else {

      
      int iarg = 1;
      
      // read commands
      bool do_stat = false;
      bool do_plot = false;
      bool do_lanczos = false;
      bool do_spectrum = false;
      bool do_output = false;
      bool do_output_mm = false;
      bool do_output_hb = false;
      bool do_output_mat = false;
      bool do_output_ghm = false;

      while (iarg<argc) {
      if (!strcmp(argv[iarg],"stat")) {
        do_stat = true;
        iarg++;
        continue;
      }
      #ifdef ESMAC_USE_PNG
      if (!strcmp(argv[iarg],"plot")) {
        do_plot = true;
        iarg++;
        continue;
      }
      #endif
      #ifdef ESMAC_USE_LAPACK
      if (!strcmp(argv[iarg],"lanczos")) {
        do_lanczos = true;
        iarg++;
        continue;
      }
      #endif
      #ifdef ESMAC_USE_LAPACK
      if (!strcmp(argv[iarg],"spectrum")) {
        do_spectrum = true;
        iarg++;
        continue;
      }
      #endif
      if (!strcmp(argv[iarg],"output")) {
        do_output = true;
        iarg++;
        continue;
      }
      if (!strcmp(argv[iarg],"mm")) {
	if (!do_output) {
	  printf("Output format specifiers [mm hb mat ghm] must be preceded by output command.\n");
	  exit(EXIT_FAILURE);
	}
        do_output_mm = true;
        iarg++;
        continue;
      }
      if (!strcmp(argv[iarg],"hb")) {
	if (!do_output) {
	  printf("Output format specifiers [mm hb mat ghm] must be preceded by output command.\n");
	  exit(EXIT_FAILURE);
	}
        do_output_hb = true;
        iarg++;
        continue;
      }
      if (!strcmp(argv[iarg],"mat")) {
	if (!do_output) {
	  printf("Output format specifiers [mm hb mat ghm] must be preceded by output command.\n");
	  exit(EXIT_FAILURE);
	}
        do_output_mat = true;
        iarg++;
        continue;
      }
      if (!strcmp(argv[iarg],"ghm")) {
	if (!do_output) {
	  printf("Output format specifiers [mm hb mat ghm] must be preceded by output command.\n");
	  exit(EXIT_FAILURE);
	}
        do_output_ghm = true;
        iarg++;
        continue;
      }
      break;
      }
      
      if (!(do_stat || do_plot || do_lanczos || do_spectrum || do_output)) {
        printf("No command given\n");
        exit(EXIT_FAILURE);
      }

      if (do_output) {
	if (!(do_output_mm || do_output_hb || do_output_mat || do_output_ghm)) {
	  printf("No output matrix format specified -> default to MatrixMarket.\n");
	  do_output_mm = true;
	}
      }

      int iarg_last_command = iarg;
      
      // read options
      bool opt_force = false;
      bool opt_quiet = false;
      bool opt_progress = false;
      bool opt_maxdim = false;
      esmac_idx_t maxdim;
      const esmac_idx_t maxdim_spectrum = 1000;
      const esmac_idx_t maxdim_lanczos = 1000000;
      bool opt_outfile = false;
      char * outfile_prefix = NULL;
      
      while (iarg<argc) {
        if (!strcmp(argv[iarg],"--force")) {
          opt_force = true;
          iarg++;
          continue;
        }
        if (!strcmp(argv[iarg],"--quiet") || !strcmp(argv[iarg],"-q")) {
          opt_quiet = true;
          iarg++;
          continue;
        }
        if (!strcmp(argv[iarg],"--progress") || !strcmp(argv[iarg],"-p")) {
          opt_progress = true;
          iarg++;
          continue;
        }
        if (!strcmp(argv[iarg],"--maxdim")) {
          opt_maxdim = true;
          iarg++;
          if (iarg >= argc) {
            printf("option \"maxdim\" requires argument\n");
            exit(EXIT_FAILURE);
          }
          maxdim = atol(argv[iarg]);
          if (maxdim <=0) {
            printf("option \"maxdim\" requires integer argument >0\n");
            exit(EXIT_FAILURE);
          }
          iarg++;
          continue;
        }
        if (!strcmp(argv[iarg],"--outfile") || !strcmp(argv[iarg],"-o")) {
          opt_outfile = true;
          iarg++;
          if (iarg >= argc) {
            printf("option \"outfile\" requires argument\n");
            exit(EXIT_FAILURE);
          }
          outfile_prefix = argv[iarg];
          iarg++;
          continue;
        }
        break;
      }

      if (iarg >= argc) {
        printf("command(s) ");
        int i;
        for (i=1;i<iarg_last_command;i++) {
          printf("\"%s\" ",argv[i]);
        }
        printf("require(s) example\n");
        exit(EXIT_FAILURE);
      }

      int filename_len;
      if (outfile_prefix) {
        filename_len = strlen(outfile_prefix);
      } else {
        filename_len = 5; // "esmac"
      }
      filename_len += 20;  // suffices
      char * filename = malloc(filename_len * sizeof *filename);

      my_gen = esmac_collection_parse(argv[iarg], &info);

      if (!opt_quiet) {
	printf("Example\n-------\n%s\n",my_gen->name);
	
	char * data_s = esmac_collection_example_data(my_gen);
      
	printf("\nParameters\n----------\n%s\n\n",data_s);
      }
      
      esmac_matrix_statistics_t st;

      esmac_idx_t nrow, ncol;
      esmac_generator_work_t * my_work = esmac_generator_alloc(my_gen, &info);
      nrow=esmac_generator_query(my_gen,my_work,"nrow");
      ncol=esmac_generator_query(my_gen,my_work,"ncol");
      esmac_generator_free(my_work);
      
      int coarse_pattern_px = 2000;
      int * coarse_pattern = NULL;

      if (do_plot) {
	coarse_pattern = malloc(coarse_pattern_px * coarse_pattern_px * sizeof *coarse_pattern );
      }

      if (do_stat || do_plot) {
	clock_t start_clock = clock();
      
	esmac_collect_matrix_statistics(my_gen, &st, coarse_pattern_px, coarse_pattern, opt_progress);
      
	clock_t diff_clock = clock() - start_clock;

	if (do_stat) {
	  if (opt_outfile) {
	    snprintf(filename,filename_len,"%s.stat.txt",outfile_prefix);
	    FILE * f = fopen(filename,"w");
	    esmac_print_statistics(f,&st);
	    fclose(f);
	  } else {
	    esmac_print_statistics(stdout, &st);
	  }
	}

  #ifdef ESMAC_USE_PNG
	if (do_plot) {
	  if (opt_outfile) {
	    snprintf(filename,filename_len,"%s.pattern.png",outfile_prefix);
	  } else {
	    snprintf(filename,filename_len,"esmacapp.pattern.png");
	  }
	  esmac_plot_pattern(coarse_pattern_px, coarse_pattern_px, coarse_pattern, 5, filename);
	}
  if (coarse_pattern) {free(coarse_pattern); coarse_pattern=NULL;}
  #endif

	if (!opt_quiet) {
	  double gen_time = (double) diff_clock / (double) CLOCKS_PER_SEC;
	  int msec = diff_clock * 1000 / CLOCKS_PER_SEC;
	  printf("\n\n====================\nGenerated %"ESMACPRIDX" rows in %d seconds %d milliseconds [%d rows/second]\n", st.nr, msec/1000, msec%1000, (int) ((double) st.nr/gen_time));  
	}
      }
      
      if (do_lanczos || do_spectrum) {
	if (nrow != ncol) {
	  printf("\"lanczos\" and \"spectrum\" require a square matrix\n");
	  exit(EXIT_FAILURE);
	}
	if (opt_maxdim) {
	  if (nrow > maxdim) {
	    printf("matrix (dim=%"ESMACPRIDX") is larger than maxdim (=%"ESMACPRIDX")\n",nrow,maxdim);
	    exit(EXIT_FAILURE);
	  } 
	} else {
	  if (do_spectrum) {
	    if (nrow > maxdim_spectrum && !opt_force) {
	      printf("matrix (dim=%"ESMACPRIDX") is larger than maxdim for \"spectrum\" (=%"ESMACPRIDX")\nUse --force to override\n",nrow,maxdim_spectrum);
	      exit(EXIT_FAILURE);
	    }
	  }
	  if (do_lanczos) {
	    if (nrow > maxdim_lanczos && !opt_force) {
	      printf("matrix (dim=%"ESMACPRIDX") is larger than maxdim for \"lanczos\" (=%"ESMACPRIDX")\nUse --force to override\n",nrow,maxdim_lanczos);
	      exit(EXIT_FAILURE);
	    }
	  }
	}
      }

      if (do_lanczos || do_spectrum || do_output) {

	esmac_sparsemat_t *sm = esmac_sparsemat_from_generator(my_gen);
  #ifdef ESMAC_USE_LAPACK
	if (do_lanczos) {
	  double tol = 1e-4;
	  double ev1,ev2, eps1,eps2;
	  info = esmac_lanczos_ev_mat(sm, tol, &ev1, &ev2, &eps1, &eps2);
	  if (opt_outfile) {
	    snprintf(filename,filename_len,"%s.lanczos.txt",outfile_prefix);
	    FILE * f = fopen(filename,"w");
	    fprintf(f,"# Lanczos for minimal & maximal eigenvalue\n");
	    fprintf(f,"%f [%f]\n",ev1,eps1);
	    fprintf(f,"%f [%f]\n",ev2,eps2);
	    fclose(f);
	  } else {
	    printf("Lanczos for minimal & maximal eigenvalue\n=======================\n");
	    printf("%f [%f]\n",ev1,eps1);
	    printf("%f [%f]\n",ev2,eps2);
	  }
	}
  #endif
  #ifdef ESMAC_USE_LAPACK
	if (do_spectrum) {
	  double *spec;
	  esmac_spectrum_mat(sm, &spec);
	  if (opt_outfile) {
	    snprintf(filename,filename_len,"%s.spectrum.txt",outfile_prefix);
	  } else {
	    snprintf(filename,filename_len,"esmacapp.spectrum.txt");
	  }
	  FILE * f = fopen(filename,"w");
	  int i;
	  for (i=0;i<nrow;i++) {
	    fprintf(f,"%d %e\n",i,spec[i]);
	  }
	  fclose(f);
	  if (spec) {free(spec);}
	}
  #endif
	if (do_output) {
	  if (do_output_mm) {
	    if (opt_outfile) {
	      snprintf(filename,filename_len,"%s.mm",outfile_prefix);
	    } else {
	      snprintf(filename,filename_len,"esmacapp.matrix.mm");
	    }
	    esmac_sparsemat_io_write_mm(sm, filename);
	  }
	  if (do_output_hb) {
	    if (opt_outfile) {
	      snprintf(filename,filename_len,"%s.hb",outfile_prefix);
	    } else {
	      snprintf(filename,filename_len,"esmacapp.matrix.hb");
	    }
	    esmac_sparsemat_io_write_hb(sm, filename);
	  }
	  if (do_output_mat) {
	    if (opt_outfile) {
	      snprintf(filename,filename_len,"%s.mat",outfile_prefix);
	    } else {
	      snprintf(filename,filename_len,"esmacapp.matrix.mat");
	    }
	    esmac_sparsemat_io_write_matlab(sm, filename);
	  }
	  if (do_output_ghm) {
	    if (opt_outfile) {
	      snprintf(filename,filename_len,"%s.ghm",outfile_prefix);
	    } else {
	      snprintf(filename,filename_len,"esmacapp.matrix.ghm");
	    }
	    esmac_sparsemat_io_write_ghost(sm, filename);
	  }
	}
	if (sm) {
	  esmac_sparsemat_free(sm);      
	}
      }
      
      esmac_collection_example_free(my_gen);

    }
    
  }
  
  return 0;
}
