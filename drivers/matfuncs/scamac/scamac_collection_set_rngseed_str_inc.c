if (!strcmp(gen->name,"Anderson")) {
  if (!strcmp(parname,"seed")) {
    ( (scamac_matrix_Anderson_params_st *) gen->par)->seed = scamac_rng_string_to_seed(seedstr);
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Exciton")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeBosonChain")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Harmonic")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Hubbard")) {
  if (!strcmp(parname,"seed")) {
    ( (scamac_matrix_Hubbard_params_st *) gen->par)->seed = scamac_rng_string_to_seed(seedstr);
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"OneFermion")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"SpinChainXXZ")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Tridiagonal")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"TridiagonalComplex")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"TridiagonalReal")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

