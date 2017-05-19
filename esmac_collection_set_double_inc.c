if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"t")) {
    ( (esmac_matrix_FreeBosonChain_params_t *) gen->par)->t = val;
    return 0;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"t")) {
    ( (esmac_matrix_FreeFermionChain_params_t *) gen->par)->t = val;
    return 0;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Harmonic")) {
  if (!strcmp(parname,"omega")) {
    ( (esmac_matrix_Harmonic_params_t *) gen->par)->omega = val;
    return 0;
  }

  if (!strcmp(parname,"lambda")) {
    ( (esmac_matrix_Harmonic_params_t *) gen->par)->lambda = val;
    return 0;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Hubbard")) {
  if (!strcmp(parname,"t")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->t = val;
    return 0;
  }

  if (!strcmp(parname,"U")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->U = val;
    return 0;
  }

  if (!strcmp(parname,"ranpot")) {
    ( (esmac_matrix_Hubbard_params_t *) gen->par)->ranpot = val;
    return 0;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"SpinChainXXZ")) {
  if (!strcmp(parname,"Jxy")) {
    ( (esmac_matrix_SpinChainXXZ_params_t *) gen->par)->Jxy = val;
    return 0;
  }

  if (!strcmp(parname,"Jz")) {
    ( (esmac_matrix_SpinChainXXZ_params_t *) gen->par)->Jz = val;
    return 0;
  }

  if (!strcmp(parname,"Bz")) {
    ( (esmac_matrix_SpinChainXXZ_params_t *) gen->par)->Bz = val;
    return 0;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

