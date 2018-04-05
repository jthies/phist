if (!strcmp(matname,"Anderson")) {
if (desc) {
 char * my_desc = malloc(76 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\nAnderson model of localization in 1D, 2D, 3D",75);
 my_desc[75]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"Exciton")) {
if (desc) {
 char * my_desc = malloc(55 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: hermitian complex\n\nExciton on a lattice",54);
 my_desc[54]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_COMPLEX;}
if (symmetry) {*symmetry=SCAMAC_HERMITIAN;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"FreeBosonChain")) {
if (desc) {
 char * my_desc = malloc(54 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\nfree bosons on a chain",53);
 my_desc[53]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"FreeFermionChain")) {
if (desc) {
 char * my_desc = malloc(56 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\nfree fermions on a chain",55);
 my_desc[55]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"Harmonic")) {
if (desc) {
 char * my_desc = malloc(59 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\nquantum Harmonic oscillator",58);
 my_desc[58]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"Hubbard")) {
if (desc) {
 char * my_desc = malloc(80 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\nHubbard models are fermionic\nsolid state models",79);
 my_desc[79]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"OneFermion")) {
if (desc) {
 char * my_desc = malloc(54 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\none fermion on a chain",53);
 my_desc[53]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"SpinChainXXZ")) {
if (desc) {
 char * my_desc = malloc(57 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: symmetric real\n\none-dimensional XXZ model",56);
 my_desc[56]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_SYMMETRIC;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"Tridiagonal")) {
if (desc) {
 char * my_desc = malloc(53 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: hermitian complex\n\nTridiagonal matrix",52);
 my_desc[52]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_COMPLEX;}
if (symmetry) {*symmetry=SCAMAC_HERMITIAN;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"TridiagonalComplex")) {
if (desc) {
 char * my_desc = malloc(65 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: general complex\n\nNon-symmetric tridiagonal matrix",64);
 my_desc[64]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_COMPLEX;}
if (symmetry) {*symmetry=SCAMAC_GENERAL;}
return SCAMAC_EOK;
}
if (!strcmp(matname,"TridiagonalReal")) {
if (desc) {
 char * my_desc = malloc(62 * sizeof *my_desc);
 strncpy(my_desc, "matrix type: general real\n\nNon-symmetric tridiagonal matrix",61);
 my_desc[61]=0;
 *desc=my_desc;
}
if (valtype) {*valtype=SCAMAC_VAL_REAL;}
if (symmetry) {*symmetry=SCAMAC_GENERAL;}
return SCAMAC_EOK;
}
