if (!strcmp(name,"FreeBosonChain")) {
  if (!strcmp(parname,"t")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_species")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_sites")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_bosons")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"boundary_conditions")) {
    return ESMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(name,"FreeFermionChain")) {
  if (!strcmp(parname,"t")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_species")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_sites")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_fermions")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"boundary_conditions")) {
    return ESMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(name,"Harmonic")) {
  if (!strcmp(parname,"omega")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"lambda")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_bos")) {
    return ESMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(name,"Hubbard")) {
  if (!strcmp(parname,"t")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"U")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_sites")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_fermions")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"boundary_conditions")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"ranpot")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"seed")) {
    return ESMAC_PAR_INT;
  }
return -1;
}

if (!strcmp(name,"SpinChainXXZ")) {
  if (!strcmp(parname,"Jxy")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"Jz")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"Bz")) {
    return ESMAC_PAR_DOUBLE;
  }
  if (!strcmp(parname,"n_sites")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"n_up")) {
    return ESMAC_PAR_INT;
  }
  if (!strcmp(parname,"boundary_conditions")) {
    return ESMAC_PAR_INT;
  }
return -1;
}

