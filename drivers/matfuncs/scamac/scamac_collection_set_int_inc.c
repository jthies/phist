if (!strcmp(gen->name,"Anderson")) {
  if (!strcmp(parname,"Lx")) {
    ( (scamac_matrix_Anderson_params_st *) gen->par)->Lx = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"Ly")) {
    ( (scamac_matrix_Anderson_params_st *) gen->par)->Ly = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"Lz")) {
    ( (scamac_matrix_Anderson_params_st *) gen->par)->Lz = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Exciton")) {
  if (!strcmp(parname,"L")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->L = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"n_species")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->n_species = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"n_sites")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->n_sites = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"n_bosons")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->n_bosons = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"n_species")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->n_species = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"n_sites")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->n_sites = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"n_fermions")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->n_fermions = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Harmonic")) {
  if (!strcmp(parname,"n_bos")) {
    ( (scamac_matrix_Harmonic_params_st *) gen->par)->n_bos = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Hubbard")) {
  if (!strcmp(parname,"n_sites")) {
    ( (scamac_matrix_Hubbard_params_st *) gen->par)->n_sites = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"n_fermions")) {
    ( (scamac_matrix_Hubbard_params_st *) gen->par)->n_fermions = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"OneFermion")) {
  if (!strcmp(parname,"n_sites")) {
    ( (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par)->n_sites = val;
//    gen->fct_unwrap_par(gen->wrapped_par, gen->par);
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"SpinChainXXZ")) {
  if (!strcmp(parname,"n_sites")) {
    ( (scamac_matrix_SpinChainXXZ_params_st *) gen->par)->n_sites = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"n_up")) {
    ( (scamac_matrix_SpinChainXXZ_params_st *) gen->par)->n_up = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"boundary_conditions")) {
    ( (scamac_matrix_SpinChainXXZ_params_st *) gen->par)->boundary_conditions = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Tridiagonal")) {
  if (!strcmp(parname,"n")) {
    ( (scamac_matrix_Tridiagonal_params_st *) gen->par)->n = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"TridiagonalComplex")) {
  if (!strcmp(parname,"n")) {
    ( (scamac_matrix_TridiagonalComplex_params_st *) gen->par)->n = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"TridiagonalReal")) {
  if (!strcmp(parname,"n")) {
    ( (scamac_matrix_TridiagonalReal_params_st *) gen->par)->n = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

