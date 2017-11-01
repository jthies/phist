if (!strcmp(matname,"FreeBosonChain")) {
  if (!strcmp(parname,"t")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_species")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_sites")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_bosons")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"PBC")) {
    return SCAMAC_PAR_BOOL;
  }
return -1;
}

if (!strcmp(matname,"FreeFermionChain")) {
  if (!strcmp(parname,"t")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_species")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_sites")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_fermions")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"PBC")) {
    return SCAMAC_PAR_BOOL;
  }
return -1;
}

if (!strcmp(matname,"Harmonic")) {
  if (!strcmp(parname,"omega")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"lambda")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_bos")) {
    return SCAMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(matname,"Hubbard")) {
  if (!strcmp(parname,"t")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"U")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_sites")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_fermions")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"boundary_conditions")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"ranpot")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"seed")) {
    return SCAMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(matname,"OneFermion")) {
  if (!strcmp(parname,"t")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_sites")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"PBC")) {
    return SCAMAC_PAR_BOOL;
  }
return -1;
}

if (!strcmp(matname,"SpinChainXXZ")) {
  if (!strcmp(parname,"Jxy")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"Jz")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"Bz")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_sites")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_up")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"boundary_conditions")) {
    return SCAMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(matname,"Tridiagonal")) {
  if (!strcmp(parname,"n")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"diag")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"offdiag")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"phi")) {
    return SCAMAC_PAR_DOUBLE;
  }
return -1;
}

if (!strcmp(matname,"TridiagonalComplex")) {
  if (!strcmp(parname,"n")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"diag_re")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"diag_im")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"supdiag")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"subdiag")) {
    return SCAMAC_PAR_DOUBLE;
  }
return -1;
}

if (!strcmp(matname,"TridiagonalReal")) {
  if (!strcmp(parname,"n")) {
    return SCAMAC_PAR_INT;
  }
  if (!strcmp(parname,"diag")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"supdiag")) {
    return SCAMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"subdiag")) {
    return SCAMAC_PAR_DOUBLE;
  }
return -1;
}

