my_string = malloc(62* sizeof *my_string);
strncpy(&my_string[0], "FreeBosonChain ",15);
strncpy(&my_string[15], "FreeFermionChain ",17);
strncpy(&my_string[32], "Harmonic ",9);
strncpy(&my_string[41], "Hubbard ",8);
strncpy(&my_string[49], "SpinChainXXZ ",13);
my_string[61] = 0;
