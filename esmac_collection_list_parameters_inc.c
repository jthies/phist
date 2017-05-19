if (!strcmp(name,"FreeBosonChain")) {
my_string = malloc(204* sizeof *my_string);
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "int n_species [number of bosonic species]\n",42);
strncpy(&my_string[70], "int n_sites [number of sites]\n",30);
strncpy(&my_string[100], "int n_bosons [number of bosons per species]\n",44);
strncpy(&my_string[144], "int boundary_conditions [ESMAC_OBC or ESMAC_PBC]\n",49);
my_string[192]=0;
return(my_string);
}

if (!strcmp(name,"FreeFermionChain")) {
my_string = malloc(210* sizeof *my_string);
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "int n_species [number of fermionic species]\n",44);
strncpy(&my_string[72], "int n_sites [number of sites]\n",30);
strncpy(&my_string[102], "int n_fermions [number of fermions per species]\n",48);
strncpy(&my_string[150], "int boundary_conditions [ESMAC_OBC or ESMAC_PBC]\n",49);
my_string[198]=0;
return(my_string);
}

if (!strcmp(name,"Harmonic")) {
my_string = malloc(105* sizeof *my_string);
strncpy(&my_string[0], "double omega [oscillator frequency]\n",36);
strncpy(&my_string[36], "double lambda [oscillator shift]\n",33);
strncpy(&my_string[69], "int n_bos [number of bosons]\n",29);
my_string[97]=0;
return(my_string);
}

if (!strcmp(name,"Hubbard")) {
my_string = malloc(292* sizeof *my_string);
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "double U [Hubbard interaction]\n",31);
strncpy(&my_string[59], "int n_sites [number of sites]\n",30);
strncpy(&my_string[89], "int n_fermions [number of fermions per spin orientation]\n",57);
strncpy(&my_string[146], "int boundary_conditions [ESMAC_OBC or ESMAC_PBC]\n",49);
strncpy(&my_string[195], "double ranpot [random on-site potential [-ranpot, ranpot]]\n",59);
strncpy(&my_string[254], "int seed [random seed]\n",23);
my_string[276]=0;
return(my_string);
}

if (!strcmp(name,"SpinChainXXZ")) {
my_string = malloc(176* sizeof *my_string);
strncpy(&my_string[0], "double Jxy [J_x=J_y]\n",21);
strncpy(&my_string[21], "double Jz [J_z]\n",16);
strncpy(&my_string[37], "double Bz [Bz]\n",15);
strncpy(&my_string[52], "int n_sites [number of sites]\n",30);
strncpy(&my_string[82], "int n_up [number of _up_ spins]\n",32);
strncpy(&my_string[114], "int boundary_conditions [ESMAC_OBC or ESMAC_PBC]\n",49);
my_string[162]=0;
return(my_string);
}

