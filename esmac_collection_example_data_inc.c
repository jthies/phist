if (!strcmp(gen->name,"FreeBosonChain")) {
  esmac_matrix_FreeBosonChain_params_t * my_FreeBosonChain_par = (esmac_matrix_FreeBosonChain_params_t *) gen->par;
 my_string = malloc(104 * sizeof *my_string);
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_FreeBosonChain_par->t);
snprintf(&my_string[strlen(my_string)],19,"%s=%d\n","n_species",my_FreeBosonChain_par->n_species);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_FreeBosonChain_par->n_sites);
snprintf(&my_string[strlen(my_string)],18,"%s=%d\n","n_bosons",my_FreeBosonChain_par->n_bosons);
snprintf(&my_string[strlen(my_string)],29,"%s=%d\n","boundary_conditions",my_FreeBosonChain_par->boundary_conditions);
return my_string;
}
if (!strcmp(gen->name,"FreeFermionChain")) {
  esmac_matrix_FreeFermionChain_params_t * my_FreeFermionChain_par = (esmac_matrix_FreeFermionChain_params_t *) gen->par;
 my_string = malloc(106 * sizeof *my_string);
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_FreeFermionChain_par->t);
snprintf(&my_string[strlen(my_string)],19,"%s=%d\n","n_species",my_FreeFermionChain_par->n_species);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_FreeFermionChain_par->n_sites);
snprintf(&my_string[strlen(my_string)],20,"%s=%d\n","n_fermions",my_FreeFermionChain_par->n_fermions);
snprintf(&my_string[strlen(my_string)],29,"%s=%d\n","boundary_conditions",my_FreeFermionChain_par->boundary_conditions);
return my_string;
}
if (!strcmp(gen->name,"Harmonic")) {
  esmac_matrix_Harmonic_params_t * my_Harmonic_par = (esmac_matrix_Harmonic_params_t *) gen->par;
 my_string = malloc(66 * sizeof *my_string);
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],25,"%s=%f\n","omega",my_Harmonic_par->omega);
snprintf(&my_string[strlen(my_string)],26,"%s=%f\n","lambda",my_Harmonic_par->lambda);
snprintf(&my_string[strlen(my_string)],15,"%s=%d\n","n_bos",my_Harmonic_par->n_bos);
return my_string;
}
if (!strcmp(gen->name,"Hubbard")) {
  esmac_matrix_Hubbard_params_t * my_Hubbard_par = (esmac_matrix_Hubbard_params_t *) gen->par;
 my_string = malloc(148 * sizeof *my_string);
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","t",my_Hubbard_par->t);
snprintf(&my_string[strlen(my_string)],21,"%s=%f\n","U",my_Hubbard_par->U);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_Hubbard_par->n_sites);
snprintf(&my_string[strlen(my_string)],20,"%s=%d\n","n_fermions",my_Hubbard_par->n_fermions);
snprintf(&my_string[strlen(my_string)],29,"%s=%d\n","boundary_conditions",my_Hubbard_par->boundary_conditions);
snprintf(&my_string[strlen(my_string)],26,"%s=%f\n","ranpot",my_Hubbard_par->ranpot);
snprintf(&my_string[strlen(my_string)],14,"%s=%d\n","seed",my_Hubbard_par->seed);
return my_string;
}
if (!strcmp(gen->name,"SpinChainXXZ")) {
  esmac_matrix_SpinChainXXZ_params_t * my_SpinChainXXZ_par = (esmac_matrix_SpinChainXXZ_params_t *) gen->par;
 my_string = malloc(127 * sizeof *my_string);
 my_string[0]=0;snprintf(&my_string[strlen(my_string)],23,"%s=%f\n","Jxy",my_SpinChainXXZ_par->Jxy);
snprintf(&my_string[strlen(my_string)],22,"%s=%f\n","Jz",my_SpinChainXXZ_par->Jz);
snprintf(&my_string[strlen(my_string)],22,"%s=%f\n","Bz",my_SpinChainXXZ_par->Bz);
snprintf(&my_string[strlen(my_string)],17,"%s=%d\n","n_sites",my_SpinChainXXZ_par->n_sites);
snprintf(&my_string[strlen(my_string)],14,"%s=%d\n","n_up",my_SpinChainXXZ_par->n_up);
snprintf(&my_string[strlen(my_string)],29,"%s=%d\n","boundary_conditions",my_SpinChainXXZ_par->boundary_conditions);
return my_string;
}
