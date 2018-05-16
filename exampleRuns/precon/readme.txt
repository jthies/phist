This folder contains some examples for preconditioner setups.
The preconditioners supported depend on the kernel library and
the TPLs found. 

Several drivers can accept a preconditioner (e.g. phist_Xbelos_blockgmres) as a pair of two strings <precon> and 
<prec_opts>. Similarly, subspacejada can read the two parameter entries "preconType" and "preconOpts" from the jadaOpts 
input file. The example parameter files are valid for both Hermitian and non-Hermitian problems, where a separate file 
for the AMG preconditioner ML for the nonsymmetric case (with the _NH_ tag in the name) is available.

This table gives an overview:

kernel lib      TPL found       precon    prec_opt

epetra          ML              ML        ml_options.xml, ml_NH_options.xml
epetra          Ifpack          ifpack    ifpack_options.xml
tpetra          Ifpack2         ifpack    ifpack_options.xml

Note that these files are meant as examples and starting point for writing your own, and not tuned to solve a particular 
problem.
