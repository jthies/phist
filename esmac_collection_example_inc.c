if (!strcmp(name,"FreeBosonChain")) {
  esmac_matrix_FreeBosonChain_params_t * my_FreeBosonChain_p = malloc(sizeof * my_FreeBosonChain_p);
  // set default values
  my_FreeBosonChain_p -> t =  1.0 ;
  my_FreeBosonChain_p -> n_species =  1 ;
  my_FreeBosonChain_p -> n_sites =  10 ;
  my_FreeBosonChain_p -> n_bosons =  5 ;
  my_FreeBosonChain_p -> boundary_conditions = 0 ;
  my_gen -> par = my_FreeBosonChain_p;
  my_gen -> gen_alloc = (esmac_gen_alloc_t) esmac_matrix_FreeBosonChain_alloc;
  my_gen -> gen_free = (esmac_gen_free_t) esmac_matrix_FreeBosonChain_free;
  my_gen -> gen_set_info = (esmac_gen_set_info_t) esmac_matrix_FreeBosonChain_set_info;
  my_gen -> gen_row = (esmac_gen_row_t) esmac_matrix_FreeBosonChain_row;
  return my_gen;
}

if (!strcmp(name,"FreeFermionChain")) {
  esmac_matrix_FreeFermionChain_params_t * my_FreeFermionChain_p = malloc(sizeof * my_FreeFermionChain_p);
  // set default values
  my_FreeFermionChain_p -> t =  1.0 ;
  my_FreeFermionChain_p -> n_species =  1 ;
  my_FreeFermionChain_p -> n_sites =  10 ;
  my_FreeFermionChain_p -> n_fermions =  5 ;
  my_FreeFermionChain_p -> boundary_conditions = 0 ;
  my_gen -> par = my_FreeFermionChain_p;
  my_gen -> gen_alloc = (esmac_gen_alloc_t) esmac_matrix_FreeFermionChain_alloc;
  my_gen -> gen_free = (esmac_gen_free_t) esmac_matrix_FreeFermionChain_free;
  my_gen -> gen_set_info = (esmac_gen_set_info_t) esmac_matrix_FreeFermionChain_set_info;
  my_gen -> gen_row = (esmac_gen_row_t) esmac_matrix_FreeFermionChain_row;
  return my_gen;
}

if (!strcmp(name,"Harmonic")) {
  esmac_matrix_Harmonic_params_t * my_Harmonic_p = malloc(sizeof * my_Harmonic_p);
  // set default values
  my_Harmonic_p -> omega =  1.0 ;
  my_Harmonic_p -> lambda =  0.0 ;
  my_Harmonic_p -> n_bos =  100 ;
  my_gen -> par = my_Harmonic_p;
  my_gen -> gen_alloc = (esmac_gen_alloc_t) esmac_matrix_Harmonic_alloc;
  my_gen -> gen_free = (esmac_gen_free_t) esmac_matrix_Harmonic_free;
  my_gen -> gen_set_info = (esmac_gen_set_info_t) esmac_matrix_Harmonic_set_info;
  my_gen -> gen_row = (esmac_gen_row_t) esmac_matrix_Harmonic_row;
  return my_gen;
}

if (!strcmp(name,"Hubbard")) {
  esmac_matrix_Hubbard_params_t * my_Hubbard_p = malloc(sizeof * my_Hubbard_p);
  // set default values
  my_Hubbard_p -> t =  1.0 ;
  my_Hubbard_p -> U =  0.0 ;
  my_Hubbard_p -> n_sites =  10 ;
  my_Hubbard_p -> n_fermions =  5 ;
  my_Hubbard_p -> boundary_conditions = 0 ;
  my_Hubbard_p -> ranpot =  0.0 ;
  my_Hubbard_p -> seed =  1 ;
  my_gen -> par = my_Hubbard_p;
  my_gen -> gen_alloc = (esmac_gen_alloc_t) esmac_matrix_Hubbard_alloc;
  my_gen -> gen_free = (esmac_gen_free_t) esmac_matrix_Hubbard_free;
  my_gen -> gen_set_info = (esmac_gen_set_info_t) esmac_matrix_Hubbard_set_info;
  my_gen -> gen_row = (esmac_gen_row_t) esmac_matrix_Hubbard_row;
  return my_gen;
}

if (!strcmp(name,"SpinChainXXZ")) {
  esmac_matrix_SpinChainXXZ_params_t * my_SpinChainXXZ_p = malloc(sizeof * my_SpinChainXXZ_p);
  // set default values
  my_SpinChainXXZ_p -> Jxy =  1.0 ;
  my_SpinChainXXZ_p -> Jz =  1.0 ;
  my_SpinChainXXZ_p -> Bz =  0.0 ;
  my_SpinChainXXZ_p -> n_sites =  10 ;
  my_SpinChainXXZ_p -> n_up =  5 ;
  my_SpinChainXXZ_p -> boundary_conditions = 0 ;
  my_gen -> par = my_SpinChainXXZ_p;
  my_gen -> gen_alloc = (esmac_gen_alloc_t) esmac_matrix_SpinChainXXZ_alloc;
  my_gen -> gen_free = (esmac_gen_free_t) esmac_matrix_SpinChainXXZ_free;
  my_gen -> gen_set_info = (esmac_gen_set_info_t) esmac_matrix_SpinChainXXZ_set_info;
  my_gen -> gen_row = (esmac_gen_row_t) esmac_matrix_SpinChainXXZ_row;
  return my_gen;
}

