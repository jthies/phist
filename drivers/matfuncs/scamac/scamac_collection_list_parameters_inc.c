if (!strcmp(matname,"Anderson")) {
my_string = malloc(377* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "int Lx [dimensions of cuboid: Lx]\n",34);
strncpy(&my_string[34], "int Ly [dimensions of cuboid: Lx]\n",34);
strncpy(&my_string[68], "int Lz [dimensions of cuboid: Lx]\n",34);
strncpy(&my_string[102], "double t [hopping strength]\n",28);
strncpy(&my_string[130], "double ranpot [random on-site potential [-ranpot, ranpot]]\n",59);
strncpy(&my_string[189], "option boundary_conditions {open,periodic} [open or periodic boundary conditions]\n",82);
strncpy(&my_string[271], "rngseed seed [random seed]\n",27);
strncpy(&my_string[298], "option sweep {simple,backforth} [mode of traversal of cuboid]\n",62);
my_string[359]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"Exciton")) {
my_string = malloc(376* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double so [spin orbit]\n",23);
strncpy(&my_string[23], "double ex [exchange]\n",21);
strncpy(&my_string[44], "double mlh [mass light hole]\n",29);
strncpy(&my_string[73], "double mhh [mass heavy hole]\n",29);
strncpy(&my_string[102], "double me [mass electron]\n",26);
strncpy(&my_string[128], "double eps [dielectric constant]\n",33);
strncpy(&my_string[161], "double lc [eff. Coulomb length]\n",32);
strncpy(&my_string[193], "double kx [momentum kx]\n",24);
strncpy(&my_string[217], "double ky [momentum ky]\n",24);
strncpy(&my_string[241], "double kz [momentum kz]\n",24);
strncpy(&my_string[265], "double a [lattice constant]\n",28);
strncpy(&my_string[293], "int L [cube length]\n",20);
strncpy(&my_string[313], "option symm {para,ortho} [symmetry]\n",36);
my_string[348]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"FreeBosonChain")) {
my_string = malloc(218* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "int n_species [number of bosonic species]\n",42);
strncpy(&my_string[70], "int n_sites [number of sites]\n",30);
strncpy(&my_string[100], "int n_bosons [number of bosons per species]\n",44);
strncpy(&my_string[144], "bool PBC [open (false) or periodic (true) boundary conditions]\n",63);
my_string[206]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"FreeFermionChain")) {
my_string = malloc(224* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "int n_species [number of fermionic species]\n",44);
strncpy(&my_string[72], "int n_sites [number of sites]\n",30);
strncpy(&my_string[102], "int n_fermions [number of fermions per species]\n",48);
strncpy(&my_string[150], "bool PBC [open (false) or periodic (true) boundary conditions]\n",63);
my_string[212]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"Harmonic")) {
my_string = malloc(105* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double omega [oscillator frequency]\n",36);
strncpy(&my_string[36], "double lambda [oscillator shift]\n",33);
strncpy(&my_string[69], "int n_bos [number of bosons]\n",29);
my_string[97]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"Hubbard")) {
my_string = malloc(329* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "double U [Hubbard interaction]\n",31);
strncpy(&my_string[59], "int n_sites [number of sites]\n",30);
strncpy(&my_string[89], "int n_fermions [number of fermions per spin orientation]\n",57);
strncpy(&my_string[146], "option boundary_conditions {open,periodic} [open or periodic boundary conditions]\n",82);
strncpy(&my_string[228], "double ranpot [random on-site potential [-ranpot, ranpot]]\n",59);
strncpy(&my_string[287], "rngseed seed [random seed]\n",27);
my_string[313]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"OneFermion")) {
my_string = malloc(128* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double t [hopping strength]\n",28);
strncpy(&my_string[28], "int n_sites [number of sites]\n",30);
strncpy(&my_string[58], "bool PBC [open (false) or periodic (true) boundary conditions]\n",63);
my_string[120]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"SpinChainXXZ")) {
my_string = malloc(178* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "double Jxy [J_x=J_y]\n",21);
strncpy(&my_string[21], "double Jz [J_z]\n",16);
strncpy(&my_string[37], "double Bz [Bz]\n",15);
strncpy(&my_string[52], "int n_sites [number of sites]\n",30);
strncpy(&my_string[82], "int n_up [number of _up_ spins]\n",32);
strncpy(&my_string[114], "int boundary_conditions [SCAMAC_OBC or SCAMAC_PBC]\n",51);
my_string[164]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"Tridiagonal")) {
my_string = malloc(147* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "int n [matrix dimension]\n",25);
strncpy(&my_string[25], "double diag [diagonal element]\n",31);
strncpy(&my_string[56], "double offdiag [off-diagonal element]\n",38);
strncpy(&my_string[94], "double phi [phase for off-diagonal element]\n",44);
my_string[137]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"TridiagonalComplex")) {
my_string = malloc(195* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "int n [matrix dimension]\n",25);
strncpy(&my_string[25], "double diag_re [diagonal element]\n",34);
strncpy(&my_string[59], "double diag_im\n",15);
strncpy(&my_string[74], "double supdiag [off-diagonal element (below diagonal)]\n",55);
strncpy(&my_string[129], "double subdiag [off-diagonal element (above diagonal)]\n",55);
my_string[183]=0;
*desc = my_string;
return SCAMAC_EOK;
}

if (!strcmp(matname,"TridiagonalReal")) {
my_string = malloc(175* sizeof *my_string);
if (!my_string) {return SCAMAC_EMALLOCFAIL;}
strncpy(&my_string[0], "int n [matrix dimension]\n",25);
strncpy(&my_string[25], "double diag [diagonal element]\n",31);
strncpy(&my_string[56], "double supdiag [off-diagonal element (below diagonal)]\n",55);
strncpy(&my_string[111], "double subdiag [off-diagonal element (above diagonal)]\n",55);
my_string[165]=0;
*desc = my_string;
return SCAMAC_EOK;
}

