if (!strcmp(gen->name,"Anderson")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_matrix_Anderson_params_st *) gen->par)->t = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"ranpot")) {
    ( (scamac_matrix_Anderson_params_st *) gen->par)->ranpot = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Exciton")) {
  if (!strcmp(parname,"so")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->so = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"ex")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->ex = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"mlh")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->mlh = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"mhh")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->mhh = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"me")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->me = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"eps")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->eps = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"lc")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->lc = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"kx")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->kx = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"ky")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->ky = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"kz")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->kz = val;
    return SCAMAC_EOK;
  }

  if (!strcmp(parname,"a")) {
    ( (scamac_matrix_Exciton_params_st *) gen->par)->a = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeBosonChain")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_matrix_FreeBosonChain_params_st *) gen->par)->t = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"FreeFermionChain")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_matrix_FreeFermionChain_params_st *) gen->par)->t = val;
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
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

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
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

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"OneFermion")) {
  if (!strcmp(parname,"t")) {
    ( (scamac_wrapper_OneFermion_params_st *) gen->wrapped_par)->t = val;
//    gen->fct_unwrap_par(gen->wrapped_par, gen->par);
    return SCAMAC_EOK;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
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

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
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

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
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

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
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

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

