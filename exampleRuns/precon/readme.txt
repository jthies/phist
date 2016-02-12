This folder contains some examples for preconditioner setups.
The preconditioners supported depend on the kernel library and
the TPLs found. 

So far the only driver that accepts a preconditioner is Xbgmres, the interface to
the Belos solver BlockGMRES. It accepts the two strings <precon> and <prec_opts>.

This table gives an overview:

kernel lib      TPL found       precon    prec_opt

epetra          Ifpack          ifpack    ifpack_options.xml
tpetra          Ifpack2         ifpack    ifpack_options.xml
