/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//======================================================================
void SUBR(kpm)(TYPE(mvec_ptr) x, TYPE(sparseMat_ptr) A, _MT_ scale, _MT_ shift, int M, _MT_* mu, int* iflag)
{
#include "phist_std_typedefs.hpp"
	PHIST_ENTER_FCN(__FUNCTION__);
	
	*iflag = 0;
	 
	TYPE(mvec_ptr) y[2];

	int n_vecs;
	PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &n_vecs, iflag), *iflag);

	_ST_* mu_tmp = new _ST_[2 * M * n_vecs]();

	y[0] = x;

	phist_const_map_ptr vector_map;
	PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A, &vector_map, iflag), *iflag);
	PHIST_CHK_IERR(SUBR(mvec_create)(&y[1], vector_map, n_vecs, iflag), *iflag);

//----------------------------------------------------------------------
	// y[1] = scale * (A-shift*I) * y[0];
	PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_add_mvec)(scale, A, -shift, y[0], 0, y[1], iflag), *iflag);

	PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(y[0], y[0], mu_tmp         , iflag), *iflag);
	PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(y[1], y[0], mu_tmp + n_vecs, iflag), *iflag);
//----------------------------------------------------------------------

	for(int m = 1; m < M; m++){
//----------------------------------------------------------------------
		PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_add_mvec)(2 * scale, A, -shift, y[m&1], -1, y[(m-1)&1], iflag), *iflag);

		PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(y[ m   &1], y[m&1], mu_tmp + (2*m  )*n_vecs, iflag), *iflag);
		PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(y[(m-1)&1], y[m&1], mu_tmp + (2*m+1)*n_vecs, iflag), *iflag);
//----------------------------------------------------------------------
	}

	for(int m = 0; m < 2 * M * n_vecs; m++) mu[m] = st::real(mu_tmp[m]);

	{	// dotProducts_2_ChebMoments
		// Straight from the essex-physics repo. Inlined, though.
		_MT_ tmp;
		for(int n = 0; n < n_vecs; n++){
			tmp = 1 / mu[n];
			mu[         n] = 1;
			mu[n_vecs + n] *= tmp;
			for(int m = 2; m <= M; m++){
				mu[(2*m-2)*n_vecs+n] = 2 * mu[(2*m-2)*n_vecs+n] * tmp - mu[         n];
				mu[(2*m-1)*n_vecs+n] = 2 * mu[(2*m-1)*n_vecs+n] * tmp - mu[n_vecs + n];
			}
		}
	}

	PHIST_CHK_IERR(SUBR(mvec_delete)(y[1], iflag), *iflag);

	delete[] mu_tmp;
}
//======================================================================
void SUBR(dos)(TYPE(sparseMat_ptr) A, _MT_ scale, _MT_ shift, int M, int R, _MT_* mu, int* iflag)
{
#include "phist_std_typedefs.hpp"
	PHIST_ENTER_FCN(__FUNCTION__);
	
	bool singlevec = *iflag & PHIST_KPM_SINGLEVEC;
	*iflag = 0;
	 
	int n_iter = 1;
	int n_vecs = R;
	if(singlevec){ n_iter = R; n_vecs = 1; }

	TYPE(mvec_ptr) x;
	phist_const_map_ptr vector_map;
	PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A, &vector_map, iflag), *iflag);
	PHIST_CHK_IERR(SUBR(mvec_create)(&x, vector_map, n_vecs, iflag), *iflag);

	_MT_* mu_tmp_r = new _MT_[2 * M * n_vecs]();

	for(int i = 0; i < n_iter; i++){
		PHIST_CHK_IERR(SUBR(mvec_random)(x, iflag), *iflag);
		PHIST_CHK_IERR(SUBR(kpm)(x, A, scale, shift, M, mu_tmp_r, iflag), *iflag);
		for(int m = 0; m < 2 * M; m++) for(int r = 0; r < n_vecs; r++) mu[m] += mu_tmp_r[m * n_vecs + r];
	}

	PHIST_CHK_IERR(SUBR(mvec_delete)(x, iflag), *iflag);

	delete[] mu_tmp_r;
}
//======================================================================
