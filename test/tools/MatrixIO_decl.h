// given a file base like xyz, tries to find a matrix file [D,S,Z,C]xyz.[mm,bin,...]
// that the kernel library can read and reads the sparse matrix.
void SUBR(read_mat)(const char* filebase,const_comm_ptr_t comm,int nglob,TYPE(sparseMat_ptr) *ptr, int* iflag);
