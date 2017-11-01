if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"PBC")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->PBC = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"PBC")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->PBC = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Harmonic")) {
  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Hubbard")) {
  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"OneFermion")) {
  if (!strcmp(parname,"PBC")) {
    ( (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par)->PBC = val;
//    gen->fct_unwrap_par(gen->wrapped_par, gen->par);
    return SCAMAC_EOK;
  }

  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"SpinChainXXZ")) {
  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Tridiagonal")) {
  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"TridiagonalComplex")) {
  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"TridiagonalReal")) {
  printf("%s: Unknown bool parameter\n",__func__);
  exit(EXIT_FAILURE);
}

