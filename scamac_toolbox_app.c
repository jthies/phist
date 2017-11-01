#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#include "scamac.h"
#include "scamac_tools.h"

const char * APPNAME = "scamact";
const char * COLLNAME = "scamac";

void print_help() {
  printf("%s - Toolbox for exploring the ScaMaC (Scalable Matrix Collection)\n\n",APPNAME);
  printf("SYNOPSIS\n  %s command(s) [options] [EXAMPLE] [PARAMETERS]\n\n",APPNAME);
  printf("Example usage:\n");
  printf("  %s list\n",APPNAME);
  printf("  %s info EXAMPLE\n",APPNAME);
  printf("  %s query EXAMPLE [PARAMETERS]\n",APPNAME);
  printf("  %s stat EXAMPLE [PARAMETERS]\n",APPNAME);
#ifdef SCAMAC_USE_PNG
  printf("  %s plot EXAMPLE [PARAMETERS]\n",APPNAME);
#endif
#ifdef SCAMAC_USE_LAPACK
  printf("  %s lanczos EXAMPLE [PARAMETERS]\n",APPNAME);
#endif
#ifdef SCAMAC_USE_LAPACK
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

  if (argc<2) {
    print_help();
  } else if (!strcmp(argv[1],"list")) {
    char * s;
    scamac_list_examples(&s);
    printf("List of examples:\n\n%s\n\n",s);
    free(s);
  } else if (!strcmp(argv[1],"info")) {
    ScamacErrorCode err = SCAMAC_EOK;
    char * sp = NULL;
    char * desc = NULL;
    if (argc<3) {
      printf("\"info\" requires matrix name. Use \"%s list\" for a list of examples.\n",APPNAME);
      err = SCAMAC_EFAIL;
    }
    if (!err) {
      err = scamac_example_parameters(argv[2], &sp);
      if (err) {
        printf("Unknown matrix example \"%s\"\nUse \"%s list\" for a list of examples.\n",argv[2],APPNAME);
      }
    }
    if (!err) {
      err = scamac_example_desc(argv[2],NULL,NULL,&desc);
      if (err) {
        printf("Unknown example \"%s\"\n",argv[2]);
      }
    }
    if (!err) {
      printf("=== Example: %s ===\n=== Description ===\n%s\n\n=== parameters ===\n",argv[2],desc);
      printf("%s\n",sp);
    }
    free(sp);
    free(desc);
  } else if (!strcmp(argv[1],"query")) {
    ScamacErrorCode err = SCAMAC_EOK;
    char * errdesc = NULL;
    ScamacGenerator * my_gen = NULL;
    if (argc<3) {
      printf("\"query\" requires example\n");
      err = SCAMAC_EFAIL;
    }
    if (!err) {
      err = scamac_parse_argstr(argv[2], &my_gen, &errdesc);
      if (err) {
        printf("Unknown example or parameters in \"%s\"\nError description: *** %s\nTry \"%s list\" or \"%s info\"\n",argv[2],errdesc,APPNAME,APPNAME);
      }
    }
    if (!err) {
      err = scamac_generator_check(my_gen,&errdesc);
      if (err) {
        printf("Problem with parameter values in \"%s\"\nError description: *** %s\n",argv[2],errdesc);
      }
    }
    if (!err) {
      SCAMAC_TRY(scamac_generator_finalize(my_gen));
      char * data_s;
      SCAMAC_TRY(scamac_generator_parameter_desc(my_gen,&data_s));
      ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
      ScamacIdx ncol = scamac_generator_query_ncol(my_gen);
      printf("=== Example: %s ===\n\n=== parameters ===\n",argv[2]);
      printf("%s\n",data_s);
      printf("dimension nrow x ncol\n%"SCAMACPRIDX" %"SCAMACPRIDX"\n",nrow,ncol);
      free(data_s);
    }
    free(errdesc);
    SCAMAC_TRY(scamac_generator_destroy(my_gen));
  } else { // commands "stat, plot, lanczos, spectrum, output"
    ScamacErrorCode err = SCAMAC_EOK;
    char * errdesc = NULL;

    ScamacGenerator * my_gen = NULL;

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

    int iarg = 1;

    while (iarg<argc) {
      if (!strcmp(argv[iarg],"stat")) {
        do_stat = true;
        iarg++;
        continue;
      }
#ifdef SCAMAC_USE_PNG
      if (!strcmp(argv[iarg],"plot")) {
        do_plot = true;
        iarg++;
        continue;
      }
#endif
#ifdef SCAMAC_USE_LAPACK
      if (!strcmp(argv[iarg],"lanczos")) {
        do_lanczos = true;
        iarg++;
        continue;
      }
#endif
#ifdef SCAMAC_USE_LAPACK
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
    int iarg_argstr = 0;

    // read options
    bool opt_force = false;
    bool opt_quiet = false;
    bool opt_progress = false;
    bool opt_maxdim = false;
    ScamacIdx maxdim = 0;
    const ScamacIdx maxdim_spectrum = 1000;
    const ScamacIdx maxdim_lanczos = 1000000;
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
      } else {
        iarg_argstr=iarg;
        iarg++;
      }
    }

    if (iarg_argstr == 0) {
      printf("command(s) ");
      int i;
      for (i=1; i<iarg_last_command; i++) {
        printf("\"%s\" ",argv[i]);
      }
      printf("require(s) example\n");
      exit(EXIT_FAILURE);
    }

    int filename_len;
    if (outfile_prefix) {
      filename_len = strlen(outfile_prefix);
    } else {
      filename_len = strlen(COLLNAME);
    }
    filename_len += 100;  // suffices
    char * filename = malloc(filename_len * sizeof *filename);

    if (!err) {
      err = scamac_parse_argstr(argv[iarg_argstr], &my_gen, &errdesc);
      if (err) {
        printf("Unknown example or parameters in \"%s\"\nError description: *** %s\nTry \"%s list\" or \"%s info\"\n",argv[iarg_argstr],errdesc,APPNAME,APPNAME);
      }
    }
    if (!err) {
      err = scamac_generator_check(my_gen,&errdesc);
      if (err) {
        printf("Problem with parameter values in \"%s\"\nError description: *** %s\n",argv[iarg_argstr],errdesc);
      }
    }

    if (!err) {
      SCAMAC_TRY(scamac_generator_finalize(my_gen));

      if (!opt_quiet) {
        printf("Example\n-------\n%s\n",scamac_generator_query_name(my_gen));

        char * data_s;
        SCAMAC_TRY(scamac_generator_parameter_desc(my_gen,&data_s));

        printf("\nParameters\n----------\n%s\n\n",data_s);
      }


      ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
      ScamacIdx ncol = scamac_generator_query_ncol(my_gen);

      if (do_stat || do_plot) {
        scamac_matrix_statistics_st st;
        scamac_matrix_pattern_st * pt = NULL;

        clock_t start_clock = clock();

        if (do_stat && !do_plot) {
          SCAMAC_TRY(scamac_collect_statistics_and_pattern(my_gen, SCAMAC_DEFAULT, &st, NULL));
        } else if (!do_stat && do_plot) {
          SCAMAC_TRY(scamac_collect_statistics_and_pattern(my_gen, SCAMAC_DEFAULT, NULL, &pt));
        } else {// do_stat && do_plot
          SCAMAC_TRY(scamac_collect_statistics_and_pattern(my_gen, SCAMAC_DEFAULT, &st, &pt));
        }

        clock_t diff_clock = clock() - start_clock;

        if (do_stat) {
          char *desc;
          scamac_statistics_print(&st, &desc);
          if (opt_outfile) {
            snprintf(filename,filename_len,"%s.stat.txt",outfile_prefix);
            FILE * f = fopen(filename,"w");
            fputs(desc,f);
            fclose(f);
          } else {
            puts(desc);
          }
          free(desc);
        }

#ifdef SCAMAC_USE_PNG
        if (do_plot) {
          /* a joke */
          char *desc;
          scamac_pattern_print(pt, &desc);
          puts(desc);
          free(desc);
          /*        */
          if (opt_outfile) {
            snprintf(filename,filename_len,"%s.pattern.png",outfile_prefix);
          } else {
            snprintf(filename,filename_len,"%s.pattern.png",COLLNAME);
          }
          scamac_plot_pattern(pt, 5, filename);
          scamac_pattern_free(pt);
        }
#endif

        if (!opt_quiet) {
          double gen_time = (double) diff_clock / (double) CLOCKS_PER_SEC;
          int msec = diff_clock * 1000 / CLOCKS_PER_SEC;
          printf("\n\n====================\nGenerated %"SCAMACPRIDX" rows in %d seconds %d milliseconds [%d rows/second]\n", nrow, msec/1000, msec%1000, (int) ((double) nrow/gen_time));
        }
      }

      if (do_lanczos || do_spectrum) {
        if (nrow != ncol) {
          printf("\"lanczos\" and \"spectrum\" require a square matrix\n");
          exit(EXIT_FAILURE);
        }
        if (opt_maxdim) {
          if (nrow > maxdim) {
            printf("matrix (dim=%"SCAMACPRIDX") is larger than maxdim (=%"SCAMACPRIDX")\n",nrow,maxdim);
            exit(EXIT_FAILURE);
          }
        } else {
          if (do_spectrum) {
            if (nrow > maxdim_spectrum && !opt_force) {
              printf("matrix (dim=%"SCAMACPRIDX") is larger than maxdim for \"spectrum\" (=%"SCAMACPRIDX")\nUse --force to override\n",nrow,maxdim_spectrum);
              exit(EXIT_FAILURE);
            }
          }
          if (do_lanczos) {
            if (nrow > maxdim_lanczos && !opt_force) {
              printf("matrix (dim=%"SCAMACPRIDX") is larger than maxdim for \"lanczos\" (=%"SCAMACPRIDX")\nUse --force to override\n",nrow,maxdim_lanczos);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

      if (do_lanczos || do_spectrum || do_output) {

        scamac_sparsemat_st *sm;
        err = scamac_sparsemat_from_generator(my_gen, &sm);
        SCAMAC_CHKERR(err);
        bool is_symmetric;
        err = scamac_sparsemat_check_symmetry(sm,NULL,&is_symmetric);
        SCAMAC_CHKERR(err);
        if (!is_symmetric) {
          printf(">>> WARNING >>> Generated matrix is not symmetric [%d].\n",is_symmetric);
        }
#ifdef SCAMAC_USE_LAPACK
        if (do_lanczos) {
          double tol = 1e-4;
          double ev1,ev2, eps1,eps2;
          err=scamac_lanczos_ev_mat(sm, tol, &ev1, &ev2, &eps1, &eps2);
          SCAMAC_CHKERR(err & ~SCAMAC_ENOTCONVERGED);
          // if (info>0) {printf("Error: Lanczos\n");exit(EXIT_FAILURE);}
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
#ifdef SCAMAC_USE_LAPACK
        if (do_spectrum) {
          double *spec;
          int flag_realspec;
          ScamacIdx valtype  = scamac_generator_query_valtype (my_gen);
          ScamacIdx symmetry = scamac_generator_query_symmetry(my_gen);
          if ( valtype == SCAMAC_VAL_REAL) {
            if (symmetry == SCAMAC_SYMMETRIC) {
              scamac_spectrum_real_symmetric(sm, &spec);
              flag_realspec=1;
            } else {
              scamac_spectrum_real_general(sm, &spec);
              flag_realspec=0;
            }
          } else {
            if (symmetry == SCAMAC_HERMITIAN) {
              scamac_spectrum_cplx_hermitian(sm, &spec);
              flag_realspec=1;
            } else {
              scamac_spectrum_cplx_general(sm, &spec);
              flag_realspec=0;
            }
          }
          if (opt_outfile) {
            snprintf(filename,filename_len,"%s.spectrum.txt",outfile_prefix);
          } else {
            snprintf(filename,filename_len,"%s.spectrum.txt",COLLNAME);
          }
          FILE * f = fopen(filename,"w");
          int i;
          if (flag_realspec) {
            for (i=0; i<nrow; i++) {
              fprintf(f,"%d %e\n",i,spec[i]);
            }
          } else {
            for (i=0; i<nrow; i++) {
              fprintf(f,"%d %e %e\n",i,spec[2*i],spec[2*i+1]);
            }
          }
          fclose(f);
          if (spec) {
            free(spec);
          }
        }
#endif
        if (do_output) {
          if (do_output_mm) {
            if (opt_outfile) {
              snprintf(filename,filename_len,"%s.mm",outfile_prefix);
            } else {
              snprintf(filename,filename_len,"%s.matrix.mm",COLLNAME);
            }
            scamac_sparsemat_io_write_mm(sm, filename);
          }
          if (do_output_hb) {
            if (opt_outfile) {
              snprintf(filename,filename_len,"%s.hb",outfile_prefix);
            } else {
              snprintf(filename,filename_len,"%s.matrix.hb",COLLNAME);
            }
            scamac_sparsemat_io_write_hb(sm, filename);
          }
          if (do_output_mat) {
            if (opt_outfile) {
              snprintf(filename,filename_len,"%s.mat",outfile_prefix);
            } else {
              snprintf(filename,filename_len,"%s.matrix.mat",COLLNAME);
            }
            scamac_sparsemat_io_write_matlab(sm, filename);
          }
          if (do_output_ghm) {
            if (opt_outfile) {
              snprintf(filename,filename_len,"%s.ghm",outfile_prefix);
            } else {
              snprintf(filename,filename_len,"%s.matrix.ghm",COLLNAME);
            }
            scamac_sparsemat_io_write_ghost(sm, filename);
          }
        }
        if (sm) {
          SCAMAC_TRY(scamac_sparsemat_free(sm));
        }
      }

    }

    SCAMAC_TRY(scamac_generator_destroy(my_gen));
    free(errdesc);

  }

  return 0;
}
