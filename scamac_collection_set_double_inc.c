if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->t = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->t = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Harmonic")) {
  if (!strcmp(parname,"omega")) {
    ( (scamac_matrix_Harmonic_params_st *) gen->par)->omega = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"lambda")) {
    ( (scamac_matrix_Harmonic_params_st *) gen->par)->lambda = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Hubbard")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_matrix_Hubbard_params_st *) gen->par)->t = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"U")) {
    ( (scamac_matrix_Hubbard_params_st *) gen->par)->U = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"ranpot")) {
    ( (scamac_matrix_Hubbard_params_st *) gen->par)->ranpot = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"OneFermion")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par)->t = val;
//    gen->fct_unwrap_par(gen->wrapped_par, gen->par);
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"SpinChainXXZ")) {
  if (!strcmp(parname,"Jxy")) {
    ( (scamac_matrix_SpinChainXXZ_params_st *) gen->par)->Jxy = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"Jz")) {
    ( (scamac_matrix_SpinChainXXZ_params_st *) gen->par)->Jz = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"Bz")) {
    ( (scamac_matrix_SpinChainXXZ_params_st *) gen->par)->Bz = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"Tridiagonal")) {
  if (!strcmp(parname,"diag")) {
    ( (scamac_matrix_Tridiagonal_params_st *) gen->par)->diag = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"offdiag")) {
    ( (scamac_matrix_Tridiagonal_params_st *) gen->par)->offdiag = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"phi")) {
    ( (scamac_matrix_Tridiagonal_params_st *) gen->par)->phi = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"TridiagonalComplex")) {
  if (!strcmp(parname,"diag_re")) {
    ( (scamac_matrix_TridiagonalComplex_params_st *) gen->par)->diag_re = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"diag_im")) {
    ( (scamac_matrix_TridiagonalComplex_params_st *) gen->par)->diag_im = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"supdiag")) {
    ( (scamac_matrix_TridiagonalComplex_params_st *) gen->par)->supdiag = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"subdiag")) {
    ( (scamac_matrix_TridiagonalComplex_params_st *) gen->par)->subdiag = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

if (!strcmp(gen->name,"TridiagonalReal")) {
  if (!strcmp(parname,"diag")) {
    ( (scamac_matrix_TridiagonalReal_params_st *) gen->par)->diag = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"supdiag")) {
    ( (scamac_matrix_TridiagonalReal_params_st *) gen->par)->supdiag = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"subdiag")) {
    ( (scamac_matrix_TridiagonalReal_params_st *) gen->par)->subdiag = val;
    return SCAMAC_EOK;
  }

  printf("%s: Unknown double parameter\n",__func__);
  exit(EXIT_FAILURE);
}

