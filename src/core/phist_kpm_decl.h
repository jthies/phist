// KPM algorithm (Weiße, Wellein, Alvermann & Fehske, Rev. Mod. Phys. 78 (1), 275-306).
// The arguments scale and shift have to map the spectrum of A into [-1,1]. 
// M is the number of iterations.
// Output array mu is assumed to have size 2 * M * #cols(x).
// Currently there are no additional options for iflag.
void SUBR(kpm)(TYPE(mvec_ptr) x, TYPE(sparseMat_ptr) A, _MT_ scale, _MT_ shift, int M, _MT_* mu, int* iflag);

// KPM-DOS algorithm (Weiße, Wellein, Alvermann & Fehske, Rev. Mod. Phys. 78 (1), 275-306).
// The arguments scale and shift have to map the spectrum of A into [-1,1]. 
// M is the number of iterations.
// R is the number of iterations.
// Output array mu is assumed to have size 2 * M.
// The flag PHIST_KPM_SINGLEVEC iterates a single vector multiple times instead of computing multiple vectors simultaneously.
void SUBR(dos)(TYPE(sparseMat_ptr) A, _MT_ scale, _MT_ shift, int M, int R, _MT_* mu, int* iflag);
