if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"n_species")) {
    ( (esmac_matrix_FreeBosonChain_params_t *) gen->par)->n_species = val;
    return 0;
  }

  if (!strcmp(parname,"n_sites")) {
    ( (esmac_matrix_FreeBosonChain_params_t *) gen->par)->n_sites = val;
    return 0;
  }

  if (!strcmp(parname,"n_bosons")) {
    ( (esmac_matrix_FreeBosonChain_params_t *) gen->par)->n_bosons = val;
    return 0;
  }

  if (!strcmp(parname,"boundary_conditions")) {
    ( (esmac_matrix_FreeBosonChain_params_t *) gen->par)->boundary_conditions = val;
    return 0;
  }

  printf("%s: Unknown integer parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"n_species")) {
    ( (esmac_matrix_FreeFermionChain_params_t *) gen->par)->n_species = val;
    return 0;
  }

  if (!strcmp(parname,"n_sites")) {
    ( (esmac_matrix_FreeFermionChain_params_t *) gen->par)->n_sites = val;
    return 0;
  }

  if (!strcmp(parname,"n_fermions")) {
    ( (esmac_matrix_FreeFermionChain_params_t *) gen->par)->n_fermions = val;
    return 0;
  }

  if (!strcmp(parname,"boundary_conditions")) {
    ( (esmac_matrix_FreeFermionChain_params_t *) gen->par)->boundary_conditions = val;
    return 0;
  }

  printf("%s: Unknown integer parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Harmonic")) {
  if (!strcmp(parname,"n_bos")) {
    ( (esmac_matrix_Harmonic_params_t *) gen->par)->n_bos = val;
    return 0;
  }

  printf("%s: Unknown integer parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Hubbard")) {
  if (!strcmp(parname,"n_sites")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->n_sites = val;
    return 0;
  }

  if (!strcmp(parname,"n_fermions")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->n_fermions = val;
    return 0;
  }

  if (!strcmp(parname,"boundary_conditions")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->boundary_conditions = val;
    return 0;
  }

  if (!strcmp(parname,"seed")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->seed = val;
    return 0;
  }

  printf("%s: Unknown integer parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"SpinChainXXZ")) {
  if (!strcmp(parname,"n_sites")) {
    ( (esmac_matrix_SpinChainXXZ_params_t *) gen->par)->n_sites = val;
    return 0;
  }

  if (!strcmp(parname,"n_up")) {
    ( (esmac_matrix_SpinChainXXZ_params_t *) gen->par)->n_up = val;
    return 0;
  }

  if (!strcmp(parname,"boundary_conditions")) {
    ( (esmac_matrix_SpinChainXXZ_params_t *) gen->par)->boundary_conditions = val;
    return 0;
  }

  printf("%s: Unknown integer parameter\n",__func__);
  exit(EXIT_FAILURE);
}

