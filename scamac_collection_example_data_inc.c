if (!strcmp(gen->name,"Anderson")) {
static const char option_Anderson_sweep[][100]={"simple","backforth"};
static const char option_Anderson_bc[][100]={"open","periodic"};
  scamac_matrix_Anderson_params_st * my_Anderson_par = (scamac_matrix_Anderson_params_st *) gen->par;
 my_string = malloc(626 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],15,format_int,"Lx",my_Anderson_par->Lx);
snprintf(&my_string[strlen(my_string)],15,format_int,"Ly",my_Anderson_par->Ly);
snprintf(&my_string[strlen(my_string)],15,format_int,"Lz",my_Anderson_par->Lz);
snprintf(&my_string[strlen(my_string)],29,format_double,"t",my_Anderson_par->t);
snprintf(&my_string[strlen(my_string)],34,format_double,"ranpot",my_Anderson_par->ranpot);
snprintf(&my_string[strlen(my_string)],122,format_option,"boundary_conditions",option_Anderson_bc[my_Anderson_par->boundary_conditions]);
snprintf(&my_string[strlen(my_string)],29,format_rngseed,"seed",my_Anderson_par->seed);
snprintf(&my_string[strlen(my_string)],108,format_option,"sweep",option_Anderson_sweep[my_Anderson_par->sweep]);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Exciton")) {
static const char option_Exciton_sym[][100]={"para","ortho"};
  scamac_matrix_Exciton_params_st * my_Exciton_par = (scamac_matrix_Exciton_params_st *) gen->par;
 my_string = malloc(712 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],30,format_double,"so",my_Exciton_par->so);
snprintf(&my_string[strlen(my_string)],30,format_double,"ex",my_Exciton_par->ex);
snprintf(&my_string[strlen(my_string)],31,format_double,"mlh",my_Exciton_par->mlh);
snprintf(&my_string[strlen(my_string)],31,format_double,"mhh",my_Exciton_par->mhh);
snprintf(&my_string[strlen(my_string)],30,format_double,"me",my_Exciton_par->me);
snprintf(&my_string[strlen(my_string)],31,format_double,"eps",my_Exciton_par->eps);
snprintf(&my_string[strlen(my_string)],30,format_double,"lc",my_Exciton_par->lc);
snprintf(&my_string[strlen(my_string)],30,format_double,"kx",my_Exciton_par->kx);
snprintf(&my_string[strlen(my_string)],30,format_double,"ky",my_Exciton_par->ky);
snprintf(&my_string[strlen(my_string)],30,format_double,"kz",my_Exciton_par->kz);
snprintf(&my_string[strlen(my_string)],29,format_double,"a",my_Exciton_par->a);
snprintf(&my_string[strlen(my_string)],14,format_int,"L",my_Exciton_par->L);
snprintf(&my_string[strlen(my_string)],107,format_option,"symm",option_Exciton_sym[my_Exciton_par->symm]);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"FreeBosonChain")) {
  scamac_matrix_FreeBosonChain_params_st * my_FreeBosonChain_par = (scamac_matrix_FreeBosonChain_params_st *) gen->par;
 my_string = malloc(367 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],29,format_double,"t",my_FreeBosonChain_par->t);
snprintf(&my_string[strlen(my_string)],22,format_int,"n_species",my_FreeBosonChain_par->n_species);
snprintf(&my_string[strlen(my_string)],20,format_int,"n_sites",my_FreeBosonChain_par->n_sites);
snprintf(&my_string[strlen(my_string)],21,format_int,"n_bosons",my_FreeBosonChain_par->n_bosons);
snprintf(&my_string[strlen(my_string)],16,format_bool,"PBC",my_FreeBosonChain_par->PBC?"True":"False");
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"FreeFermionChain")) {
  scamac_matrix_FreeFermionChain_params_st * my_FreeFermionChain_par = (scamac_matrix_FreeFermionChain_params_st *) gen->par;
 my_string = malloc(369 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],29,format_double,"t",my_FreeFermionChain_par->t);
snprintf(&my_string[strlen(my_string)],22,format_int,"n_species",my_FreeFermionChain_par->n_species);
snprintf(&my_string[strlen(my_string)],20,format_int,"n_sites",my_FreeFermionChain_par->n_sites);
snprintf(&my_string[strlen(my_string)],23,format_int,"n_fermions",my_FreeFermionChain_par->n_fermions);
snprintf(&my_string[strlen(my_string)],16,format_bool,"PBC",my_FreeFermionChain_par->PBC?"True":"False");
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Harmonic")) {
  scamac_matrix_Harmonic_params_st * my_Harmonic_par = (scamac_matrix_Harmonic_params_st *) gen->par;
 my_string = malloc(344 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],33,format_double,"omega",my_Harmonic_par->omega);
snprintf(&my_string[strlen(my_string)],34,format_double,"lambda",my_Harmonic_par->lambda);
snprintf(&my_string[strlen(my_string)],18,format_int,"n_bos",my_Harmonic_par->n_bos);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Hubbard")) {
static const char option_Hubbard_bc[][100]={"open","periodic"};
  scamac_matrix_Hubbard_params_st * my_Hubbard_par = (scamac_matrix_Hubbard_params_st *) gen->par;
 my_string = malloc(545 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],29,format_double,"t",my_Hubbard_par->t);
snprintf(&my_string[strlen(my_string)],29,format_double,"U",my_Hubbard_par->U);
snprintf(&my_string[strlen(my_string)],20,format_int,"n_sites",my_Hubbard_par->n_sites);
snprintf(&my_string[strlen(my_string)],23,format_int,"n_fermions",my_Hubbard_par->n_fermions);
snprintf(&my_string[strlen(my_string)],122,format_option,"boundary_conditions",option_Hubbard_bc[my_Hubbard_par->boundary_conditions]);
snprintf(&my_string[strlen(my_string)],34,format_double,"ranpot",my_Hubbard_par->ranpot);
snprintf(&my_string[strlen(my_string)],29,format_rngseed,"seed",my_Hubbard_par->seed);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"OneFermion")) {
  scamac_wrapper_OneFermion_params_st * my_OneFermion_par = (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par;
 my_string = malloc(324 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],29,format_double,"t",my_OneFermion_par->t);
snprintf(&my_string[strlen(my_string)],20,format_int,"n_sites",my_OneFermion_par->n_sites);
snprintf(&my_string[strlen(my_string)],16,format_bool,"PBC",my_OneFermion_par->PBC?"True":"False");
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"SpinChainXXZ")) {
  scamac_matrix_SpinChainXXZ_params_st * my_SpinChainXXZ_par = (scamac_matrix_SpinChainXXZ_params_st *) gen->par;
 my_string = malloc(419 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],31,format_double,"Jxy",my_SpinChainXXZ_par->Jxy);
snprintf(&my_string[strlen(my_string)],30,format_double,"Jz",my_SpinChainXXZ_par->Jz);
snprintf(&my_string[strlen(my_string)],30,format_double,"Bz",my_SpinChainXXZ_par->Bz);
snprintf(&my_string[strlen(my_string)],20,format_int,"n_sites",my_SpinChainXXZ_par->n_sites);
snprintf(&my_string[strlen(my_string)],17,format_int,"n_up",my_SpinChainXXZ_par->n_up);
snprintf(&my_string[strlen(my_string)],32,format_int,"boundary_conditions",my_SpinChainXXZ_par->boundary_conditions);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Tridiagonal")) {
  scamac_matrix_Tridiagonal_params_st * my_Tridiagonal_par = (scamac_matrix_Tridiagonal_params_st *) gen->par;
 my_string = malloc(371 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],14,format_int,"n",my_Tridiagonal_par->n);
snprintf(&my_string[strlen(my_string)],32,format_double,"diag",my_Tridiagonal_par->diag);
snprintf(&my_string[strlen(my_string)],35,format_double,"offdiag",my_Tridiagonal_par->offdiag);
snprintf(&my_string[strlen(my_string)],31,format_double,"phi",my_Tridiagonal_par->phi);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"TridiagonalComplex")) {
  scamac_matrix_TridiagonalComplex_params_st * my_TridiagonalComplex_par = (scamac_matrix_TridiagonalComplex_params_st *) gen->par;
 my_string = malloc(413 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],14,format_int,"n",my_TridiagonalComplex_par->n);
snprintf(&my_string[strlen(my_string)],35,format_double,"diag_re",my_TridiagonalComplex_par->diag_re);
snprintf(&my_string[strlen(my_string)],35,format_double,"diag_im",my_TridiagonalComplex_par->diag_im);
snprintf(&my_string[strlen(my_string)],35,format_double,"supdiag",my_TridiagonalComplex_par->supdiag);
snprintf(&my_string[strlen(my_string)],35,format_double,"subdiag",my_TridiagonalComplex_par->subdiag);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"TridiagonalReal")) {
  scamac_matrix_TridiagonalReal_params_st * my_TridiagonalReal_par = (scamac_matrix_TridiagonalReal_params_st *) gen->par;
 my_string = malloc(375 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;if (print_name) {
  snprintf(&my_string[strlen(my_string)],258, format_name, gen->name);
}
snprintf(&my_string[strlen(my_string)],14,format_int,"n",my_TridiagonalReal_par->n);
snprintf(&my_string[strlen(my_string)],32,format_double,"diag",my_TridiagonalReal_par->diag);
snprintf(&my_string[strlen(my_string)],35,format_double,"supdiag",my_TridiagonalReal_par->supdiag);
snprintf(&my_string[strlen(my_string)],35,format_double,"subdiag",my_TridiagonalReal_par->subdiag);
 if (print_as_argstr) { delete_last_char(sepchar, my_string); } 
*desc = my_string;
return SCAMAC_EOK;
}
