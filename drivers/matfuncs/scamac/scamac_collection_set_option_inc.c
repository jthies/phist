if (!strcmp(gen->name,"Anderson")) {
  if (!strcmp(parname,"boundary_conditions")) {
    if (!strcmp(option,"open")) {
      ( (scamac_matrix_Anderson_params_st *) gen->par)->boundary_conditions = Anderson_open;
      return SCAMAC_EOK;
    }
    if (!strcmp(option,"periodic")) {
      ( (scamac_matrix_Anderson_params_st *) gen->par)->boundary_conditions = Anderson_periodic;
      return SCAMAC_EOK;
    }
    return SCAMAC_EINVALID | 3 << SCAMAC_ESHIFT;
  }

  if (!strcmp(parname,"sweep")) {
    if (!strcmp(option,"simple")) {
      ( (scamac_matrix_Anderson_params_st *) gen->par)->sweep = Anderson_simple;
      return SCAMAC_EOK;
    }
    if (!strcmp(option,"backforth")) {
      ( (scamac_matrix_Anderson_params_st *) gen->par)->sweep = Anderson_backforth;
      return SCAMAC_EOK;
    }
    return SCAMAC_EINVALID | 3 << SCAMAC_ESHIFT;
  }

  // Unknown parameter
  return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
}

if (!strcmp(gen->name,"Exciton")) {
  if (!strcmp(parname,"symm")) {
    if (!strcmp(option,"para")) {
      ( (scamac_matrix_Exciton_params_st *) gen->par)->symm = Exciton_para;
      return SCAMAC_EOK;
    }
    if (!strcmp(option,"ortho")) {
      ( (scamac_matrix_Exciton_params_st *) gen->par)->symm = Exciton_ortho;
      return SCAMAC_EOK;
    }
    return SCAMAC_EINVALID | 3 << SCAMAC_ESHIFT;
  }

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
  if (!strcmp(parname,"boundary_conditions")) {
    if (!strcmp(option,"open")) {
      ( (scamac_matrix_Hubbard_params_st *) gen->par)->boundary_conditions = Hubbard_open;
      return SCAMAC_EOK;
    }
    if (!strcmp(option,"periodic")) {
      ( (scamac_matrix_Hubbard_params_st *) gen->par)->boundary_conditions = Hubbard_periodic;
      return SCAMAC_EOK;
    }
    return SCAMAC_EINVALID | 3 << SCAMAC_ESHIFT;
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

