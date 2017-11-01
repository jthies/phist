if (!strcmp(gen->name,"FreeBosonChain")) {
  scamac_matrix_FreeBosonChain_params_st * my_FreeBosonChain_par = (scamac_matrix_FreeBosonChain_params_st *) gen->par;
 my_string = malloc(88 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_FreeBosonChain_par->t);
snprintf(&my_string[strlen(my_string)],19,"%s=%d\n","n_species",my_FreeBosonChain_par->n_species);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_FreeBosonChain_par->n_sites);
snprintf(&my_string[strlen(my_string)],18,"%s=%d\n","n_bosons",my_FreeBosonChain_par->n_bosons);
snprintf(&my_string[strlen(my_string)],13,"%s=%s\n","PBC",my_FreeBosonChain_par->PBC?"True":"False");
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"FreeFermionChain")) {
  scamac_matrix_FreeFermionChain_params_st * my_FreeFermionChain_par = (scamac_matrix_FreeFermionChain_params_st *) gen->par;
 my_string = malloc(90 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_FreeFermionChain_par->t);
snprintf(&my_string[strlen(my_string)],19,"%s=%d\n","n_species",my_FreeFermionChain_par->n_species);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_FreeFermionChain_par->n_sites);
snprintf(&my_string[strlen(my_string)],20,"%s=%d\n","n_fermions",my_FreeFermionChain_par->n_fermions);
snprintf(&my_string[strlen(my_string)],13,"%s=%s\n","PBC",my_FreeFermionChain_par->PBC?"True":"False");
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Harmonic")) {
  scamac_matrix_Harmonic_params_st * my_Harmonic_par = (scamac_matrix_Harmonic_params_st *) gen->par;
 my_string = malloc(66 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],25,"%s=%f\n","omega",my_Harmonic_par->omega);
snprintf(&my_string[strlen(my_string)],26,"%s=%f\n","lambda",my_Harmonic_par->lambda);
snprintf(&my_string[strlen(my_string)],15,"%s=%d\n","n_bos",my_Harmonic_par->n_bos);
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Hubbard")) {
  scamac_matrix_Hubbard_params_st * my_Hubbard_par = (scamac_matrix_Hubbard_params_st *) gen->par;
 my_string = malloc(148 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_Hubbard_par->t);
snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","U",my_Hubbard_par->U);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_Hubbard_par->n_sites);
snprintf(&my_string[strlen(my_string)],20,"%s=%d\n","n_fermions",my_Hubbard_par->n_fermions);
snprintf(&my_string[strlen(my_string)],29,"%s=%d\n","boundary_conditions",my_Hubbard_par->boundary_conditions);
snprintf(&my_string[strlen(my_string)],26,"%s=%f\n","ranpot",my_Hubbard_par->ranpot);
snprintf(&my_string[strlen(my_string)],14,"%s=%d\n","seed",my_Hubbard_par->seed);
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"OneFermion")) {
  scamac_wrapper_OneFermion_params_st * my_OneFermion_par = (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par;
 my_string = malloc(51 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_OneFermion_par->t);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_OneFermion_par->n_sites);
snprintf(&my_string[strlen(my_string)],13,"%s=%s\n","PBC",my_OneFermion_par->PBC?"True":"False");
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"SpinChainXXZ")) {
  scamac_matrix_SpinChainXXZ_params_st * my_SpinChainXXZ_par = (scamac_matrix_SpinChainXXZ_params_st *) gen->par;
 my_string = malloc(127 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],23,"%s=%f\n","Jxy",my_SpinChainXXZ_par->Jxy);
snprintf(&my_string[strlen(my_string)],22,"%s=%f\n","Jz",my_SpinChainXXZ_par->Jz);
snprintf(&my_string[strlen(my_string)],22,"%s=%f\n","Bz",my_SpinChainXXZ_par->Bz);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_SpinChainXXZ_par->n_sites);
snprintf(&my_string[strlen(my_string)],14,"%s=%d\n","n_up",my_SpinChainXXZ_par->n_up);
snprintf(&my_string[strlen(my_string)],29,"%s=%d\n","boundary_conditions",my_SpinChainXXZ_par->boundary_conditions);
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"Tridiagonal")) {
  scamac_matrix_Tridiagonal_params_st * my_Tridiagonal_par = (scamac_matrix_Tridiagonal_params_st *) gen->par;
 my_string = malloc(85 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],11,"%s=%d\n","n",my_Tridiagonal_par->n);
snprintf(&my_string[strlen(my_string)],24,"%s=%f\n","diag",my_Tridiagonal_par->diag);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","offdiag",my_Tridiagonal_par->offdiag);
snprintf(&my_string[strlen(my_string)],23,"%s=%f\n","phi",my_Tridiagonal_par->phi);
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"TridiagonalComplex")) {
  scamac_matrix_TridiagonalComplex_params_st * my_TridiagonalComplex_par = (scamac_matrix_TridiagonalComplex_params_st *) gen->par;
 my_string = malloc(119 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],11,"%s=%d\n","n",my_TridiagonalComplex_par->n);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","diag_re",my_TridiagonalComplex_par->diag_re);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","diag_im",my_TridiagonalComplex_par->diag_im);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","supdiag",my_TridiagonalComplex_par->supdiag);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","subdiag",my_TridiagonalComplex_par->subdiag);
*desc = my_string;
return SCAMAC_EOK;
}
if (!strcmp(gen->name,"TridiagonalReal")) {
  scamac_matrix_TridiagonalReal_params_st * my_TridiagonalReal_par = (scamac_matrix_TridiagonalReal_params_st *) gen->par;
 my_string = malloc(89 * sizeof *my_string);
 if (!my_string) { return SCAMAC_EMALLOCFAIL; }
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],11,"%s=%d\n","n",my_TridiagonalReal_par->n);
snprintf(&my_string[strlen(my_string)],24,"%s=%f\n","diag",my_TridiagonalReal_par->diag);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","supdiag",my_TridiagonalReal_par->supdiag);
snprintf(&my_string[strlen(my_string)],27,"%s=%f\n","subdiag",my_TridiagonalReal_par->subdiag);
*desc = my_string;
return SCAMAC_EOK;
}
