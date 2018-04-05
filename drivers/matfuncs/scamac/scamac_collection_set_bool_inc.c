if (!strcmp(gen->name,"Anderson")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Exciton")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"PBC")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->PBC = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"PBC")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->PBC = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Harmonic")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Hubbard")) {
  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"OneFermion")) {
  if (!strcmp(parname,"PBC")) {
    ( (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par)->PBC = val;
//    gen->fct_unwrap_par(gen->wrapped_par, gen->par);
    return SCAMAC_EOK;
  }

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

