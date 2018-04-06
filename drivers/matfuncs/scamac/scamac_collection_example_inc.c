if (!strcmp(matname,"Anderson")) {
  scamac_matrix_Anderson_params_st * my_Anderson_p = malloc(sizeof * my_Anderson_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_Anderson_p;
  my_gen -> fct_work_alloc = (scamac_gen_work_alloc_ft) scamac_matrix_Anderson_work_alloc;
  my_gen -> fct_work_free = (scamac_gen_work_free_ft) scamac_matrix_Anderson_work_free;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_Anderson_tables_create;
  my_gen -> fct_tables_destroy = NULL;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_Anderson_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_Anderson_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_SYMMETRIC;
  // set default values
  my_Anderson_p -> Lx =  5 ;
  my_Anderson_p -> Ly =  5 ;
  my_Anderson_p -> Lz =  5 ;
  my_Anderson_p -> t =  1.0 ;
  my_Anderson_p -> ranpot =  0.0 ;
  my_Anderson_p -> boundary_conditions = 0 ;
  my_Anderson_p -> seed =  1 ;
  my_Anderson_p -> sweep =  Anderson_simple ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"Exciton")) {
  scamac_matrix_Exciton_params_st * my_Exciton_p = malloc(sizeof * my_Exciton_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_Exciton_p;
  my_gen -> fct_work_alloc = NULL;
  my_gen -> fct_work_free = NULL;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_Exciton_tables_create;
  my_gen -> fct_tables_destroy = (scamac_gen_tables_destroy_ft) scamac_matrix_Exciton_tables_destroy;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_Exciton_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_Exciton_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_COMPLEX;
  my_gen -> info.symmetry = SCAMAC_HERMITIAN;
  // set default values
  my_Exciton_p -> so =  128.0 ;
  my_Exciton_p -> ex =  666.0 ;
  my_Exciton_p -> mlh =  0.16 ;
  my_Exciton_p -> mhh =  3.1 ;
  my_Exciton_p -> me =  0.99 ;
  my_Exciton_p -> eps =  6.94 ;
  my_Exciton_p -> lc =  1.75 ;
  my_Exciton_p -> kx =  0.0 ;
  my_Exciton_p -> ky =  0.0 ;
  my_Exciton_p -> kz =  0.0 ;
  my_Exciton_p -> a =  0.42696 ;
  my_Exciton_p -> L =  10 ;
  my_Exciton_p -> symm =  Exciton_para ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"FreeBosonChain")) {
  scamac_matrix_FreeBosonChain_params_st * my_FreeBosonChain_p = malloc(sizeof * my_FreeBosonChain_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_FreeBosonChain_p;
  my_gen -> fct_work_alloc = (scamac_gen_work_alloc_ft) scamac_matrix_FreeBosonChain_work_alloc;
  my_gen -> fct_work_free = (scamac_gen_work_free_ft) scamac_matrix_FreeBosonChain_work_free;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_FreeBosonChain_tables_create;
  my_gen -> fct_tables_destroy = (scamac_gen_tables_destroy_ft) scamac_matrix_FreeBosonChain_tables_destroy;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_FreeBosonChain_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_FreeBosonChain_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_SYMMETRIC;
  // set default values
  my_FreeBosonChain_p -> t =  1.0 ;
  my_FreeBosonChain_p -> n_species =  1 ;
  my_FreeBosonChain_p -> n_sites =  10 ;
  my_FreeBosonChain_p -> n_bosons =  5 ;
  my_FreeBosonChain_p -> PBC =  true ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"FreeFermionChain")) {
  scamac_matrix_FreeFermionChain_params_st * my_FreeFermionChain_p = malloc(sizeof * my_FreeFermionChain_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_FreeFermionChain_p;
  my_gen -> fct_work_alloc = (scamac_gen_work_alloc_ft) scamac_matrix_FreeFermionChain_work_alloc;
  my_gen -> fct_work_free = (scamac_gen_work_free_ft) scamac_matrix_FreeFermionChain_work_free;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_FreeFermionChain_tables_create;
  my_gen -> fct_tables_destroy = (scamac_gen_tables_destroy_ft) scamac_matrix_FreeFermionChain_tables_destroy;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_FreeFermionChain_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_FreeFermionChain_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_SYMMETRIC;
  // set default values
  my_FreeFermionChain_p -> t =  1.0 ;
  my_FreeFermionChain_p -> n_species =  1 ;
  my_FreeFermionChain_p -> n_sites =  10 ;
  my_FreeFermionChain_p -> n_fermions =  5 ;
  my_FreeFermionChain_p -> PBC =  true ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"Harmonic")) {
  scamac_matrix_Harmonic_params_st * my_Harmonic_p = malloc(sizeof * my_Harmonic_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_Harmonic_p;
  my_gen -> fct_work_alloc = NULL;
  my_gen -> fct_work_free = NULL;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_Harmonic_tables_create;
  my_gen -> fct_tables_destroy = NULL;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_Harmonic_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_Harmonic_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_SYMMETRIC;
  // set default values
  my_Harmonic_p -> omega =  1.0 ;
  my_Harmonic_p -> lambda =  0.0 ;
  my_Harmonic_p -> n_bos =  100 ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"Hubbard")) {
  scamac_matrix_Hubbard_params_st * my_Hubbard_p = malloc(sizeof * my_Hubbard_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_Hubbard_p;
  my_gen -> fct_work_alloc = (scamac_gen_work_alloc_ft) scamac_matrix_Hubbard_work_alloc;
  my_gen -> fct_work_free = (scamac_gen_work_free_ft) scamac_matrix_Hubbard_work_free;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_Hubbard_tables_create;
  my_gen -> fct_tables_destroy = (scamac_gen_tables_destroy_ft) scamac_matrix_Hubbard_tables_destroy;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_Hubbard_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_Hubbard_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_SYMMETRIC;
  // set default values
  my_Hubbard_p -> t =  1.0 ;
  my_Hubbard_p -> U =  0.0 ;
  my_Hubbard_p -> n_sites =  10 ;
  my_Hubbard_p -> n_fermions =  5 ;
  my_Hubbard_p -> boundary_conditions = 0 ;
  my_Hubbard_p -> ranpot =  0.0 ;
  my_Hubbard_p -> seed =  1 ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"OneFermion")) {
  scamac_wrapper_OneFermion_params_st * my_OneFermion_p = malloc(sizeof * my_OneFermion_p);
  ScamacErrorCode err;
  ScamacGenerator * my_gen;
  err = scamac_generator_obtain("FreeFermionChain", &my_gen);
  if (err) {return err;}
  strncpy(my_gen->name,"OneFermion", SCAMAC_NAME_LENGTH);
  my_gen -> wrapped_par = my_OneFermion_p;
  my_gen -> fct_unwrap_par = (scamac_gen_unwrap_ft) scamac_wrapper_OneFermion_unwrap;
  my_gen -> fct_check_wrapped = (scamac_gen_check_ft) scamac_wrapper_OneFermion_check;
  my_gen -> needs_unwrapping = true;
  // set default values
  my_OneFermion_p -> t =  1.0 ;
  my_OneFermion_p -> n_sites =  10 ;
  my_OneFermion_p -> PBC =  true ;
//  my_gen -> fct_unwrap_par(my_gen->wrapped_par, my_gen->par);
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"SpinChainXXZ")) {
  scamac_matrix_SpinChainXXZ_params_st * my_SpinChainXXZ_p = malloc(sizeof * my_SpinChainXXZ_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_SpinChainXXZ_p;
  my_gen -> fct_work_alloc = (scamac_gen_work_alloc_ft) scamac_matrix_SpinChainXXZ_work_alloc;
  my_gen -> fct_work_free = (scamac_gen_work_free_ft) scamac_matrix_SpinChainXXZ_work_free;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_SpinChainXXZ_tables_create;
  my_gen -> fct_tables_destroy = (scamac_gen_tables_destroy_ft) scamac_matrix_SpinChainXXZ_tables_destroy;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_SpinChainXXZ_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_SpinChainXXZ_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_SYMMETRIC;
  // set default values
  my_SpinChainXXZ_p -> Jxy =  1.0 ;
  my_SpinChainXXZ_p -> Jz =  1.0 ;
  my_SpinChainXXZ_p -> Bz =  0.0 ;
  my_SpinChainXXZ_p -> n_sites =  10 ;
  my_SpinChainXXZ_p -> n_up =  5 ;
  my_SpinChainXXZ_p -> boundary_conditions = 0 ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"Tridiagonal")) {
  scamac_matrix_Tridiagonal_params_st * my_Tridiagonal_p = malloc(sizeof * my_Tridiagonal_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_Tridiagonal_p;
  my_gen -> fct_work_alloc = NULL;
  my_gen -> fct_work_free = NULL;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_Tridiagonal_tables_create;
  my_gen -> fct_tables_destroy = NULL;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_Tridiagonal_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_Tridiagonal_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_COMPLEX;
  my_gen -> info.symmetry = SCAMAC_HERMITIAN;
  // set default values
  my_Tridiagonal_p -> n =  100 ;
  my_Tridiagonal_p -> diag =  0.0 ;
  my_Tridiagonal_p -> offdiag =  1.0 ;
  my_Tridiagonal_p -> phi =  0.0 ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"TridiagonalComplex")) {
  scamac_matrix_TridiagonalComplex_params_st * my_TridiagonalComplex_p = malloc(sizeof * my_TridiagonalComplex_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_TridiagonalComplex_p;
  my_gen -> fct_work_alloc = NULL;
  my_gen -> fct_work_free = NULL;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_TridiagonalComplex_tables_create;
  my_gen -> fct_tables_destroy = NULL;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_TridiagonalComplex_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_TridiagonalComplex_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_COMPLEX;
  my_gen -> info.symmetry = SCAMAC_GENERAL;
  // set default values
  my_TridiagonalComplex_p -> n =  100 ;
  my_TridiagonalComplex_p -> diag_re =  0.0 ;
  my_TridiagonalComplex_p -> diag_im =  0.0 ;
  my_TridiagonalComplex_p -> supdiag =  1.0 ;
  my_TridiagonalComplex_p -> subdiag =  1.0 ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

if (!strcmp(matname,"TridiagonalReal")) {
  scamac_matrix_TridiagonalReal_params_st * my_TridiagonalReal_p = malloc(sizeof * my_TridiagonalReal_p);
  ScamacGenerator * my_gen = malloc(sizeof * my_gen);
  strncpy(my_gen->name, matname, SCAMAC_NAME_LENGTH);
  my_gen -> name[SCAMAC_NAME_LENGTH-1]=0;
  my_gen -> par = my_TridiagonalReal_p;
  my_gen -> fct_work_alloc = NULL;
  my_gen -> fct_work_free = NULL;
  my_gen -> fct_tables_create = (scamac_gen_tables_create_ft) scamac_matrix_TridiagonalReal_tables_create;
  my_gen -> fct_tables_destroy = NULL;
  my_gen -> fct_gen_row = (scamac_gen_row_ft) scamac_matrix_TridiagonalReal_generate_row;
  my_gen -> fct_check = (scamac_gen_check_ft) scamac_matrix_TridiagonalReal_check;
  my_gen -> wrapped_par = NULL;
  my_gen -> fct_unwrap_par = NULL;
  my_gen -> needs_finalization = true;
  my_gen -> needs_unwrapping = false;
  my_gen -> tables = NULL;
  my_gen -> info.nrow = 0;
  my_gen -> info.ncol = 0;
  my_gen -> info.maxnzrow = 0;
  my_gen -> info.maxnzcol = 0;
  my_gen -> info.maxnz = 0;
  my_gen -> info.valtype = SCAMAC_VAL_REAL;
  my_gen -> info.symmetry = SCAMAC_GENERAL;
  // set default values
  my_TridiagonalReal_p -> n =  100 ;
  my_TridiagonalReal_p -> diag =  0.0 ;
  my_TridiagonalReal_p -> supdiag =  1.0 ;
  my_TridiagonalReal_p -> subdiag =  1.0 ;
  *gen = my_gen;
  return SCAMAC_EOK;
}

