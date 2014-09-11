#include "matfuncs.h"

#ifdef GHOST_HAVE_LONGIDX
  #define idx_t uint64_t
#else
  #define idx_t uint32_t
#endif


int cols_compar(  const void *a, const void *b ){ 
	return (int)(*(ghost_gidx_t *)a - *(ghost_gidx_t *)b );
	}

void cols_qsort(ghost_gidx_t *cols, char * vals ,  size_t size_v ,  ghost_gidx_t n ){
	
	size_t size_i = sizeof(ghost_gidx_t);
	size_t size   = size_v + size_i;
	char * base   = malloc( n*size );
	ghost_gidx_t i;
	for(i=0;i<n;i++) {
		memcpy( base + i*size,          cols + i       , size_i);
		memcpy( base + i*size + size_i, vals + i*size_v, size_v);}
	
	qsort ( base, n, size , cols_compar);
	
	for(i=0;i<n;i++) {
		memcpy( cols + i       , base + i*size         , size_i);
		memcpy( vals + i*size_v, base + i*size + size_i, size_v);}
	
	free(base);
}

int cols_error_check(ghost_gidx_t dim, ghost_gidx_t row, ghost_gidx_t nnz,  ghost_gidx_t * cols ){
	
	ghost_gidx_t i;
	for(i=0; i<nnz; i++) if( (cols[i] >= dim) || (cols[i] < 0) ) {printf("error: row %ld: col[%ld] = %ld >= %ld \n", row, i, cols[i], dim); exit(0);}
	
	return 0;
	
	}


//#define DEBUG

#include <limits.h>

char * int2bin(int i)
{
    size_t bits = sizeof(int) * CHAR_BIT;

    char * str = malloc(bits + 1);
    if(!str) return NULL;
    str[bits] = 0;

    // type punning because signed shift is implementation-defined
    unsigned u = *(unsigned *)&i;
    for(; bits--; u >>= 1)
    	str[bits] = u & 1 ? '1' : '0';

    return str;
}

#define BIT_OR  |
#define BIT_AND &
#define BIT_XOR ^
#define BIT_NOT ~
#define BIT_SHIFT_L <<
#define BIT_SHIFT_R >>

#define bitor( a, b ) (((idx_t)(a))|((idx_t)(b)))
#define bitand( a, b ) (((idx_t)(a))&((idx_t)(b)))
#define bitxor( a, b ) (((idx_t)(a))^((idx_t)(b)))
#define bitnot( a ) (~((idx_t)(a)))
#define bitshift_l( a, b ) (((idx_t)(a))<<((idx_t)(b)))
#define bitshift_r( a, b ) (((idx_t)(a))>>((idx_t)(b)))
#define bitor( a, b ) (((idx_t)(a))|((idx_t)(b)))


/*
idx_t bitor( idx_t a,  idx_t b ){
	return a|b;
}

idx_t bitand( idx_t a,  idx_t b ){
	return a&b;
}

idx_t bitxor( idx_t a,  idx_t b ){
	return a^b;
}

idx_t bitnot( idx_t a ){
	return ~a;
}


idx_t bitshift_l( idx_t z,  int32_t n ){
	ghost_gidx_t i;
	
#ifdef GHOST_HAVE_LONGIDX_
	ghost_gidx_t r=0;
	uint32_t tmp;
	uint32_t * t0 = (uint32_t *)( &z);
	uint32_t * t1 = (uint32_t *)( &r);
	
	if (n < 32 ){
		t1[0] = t0[0] << n; t1[1] = t0[1] << n;
		tmp =  t0[0] >> (32-n);
		t1[1] |=  tmp;
	}
	else{
		tmp =  t0[0] << (n-32);
		t1[1] |=  tmp;
	}
	return r;
	
	//if (n > 0) for( i=0; i<n; i++ ) z *= 2;
	//else       for( i=0; i<n; i++ ) z /= 2;
	//return z;
#else
	return z << n;
#endif
	}

idx_t bitshift_r( idx_t z,  int32_t  n ){
	idx_t i;

#ifdef GHOST_HAVE_LONGIDX_
		idx_t r=0;
	uint32_t tmp;
	uint32_t * t0 = (uint32_t *)( &z);
	uint32_t * t1 = (uint32_t *)( &r);
	
	if (n < 32 ){
		t1[0] = t0[0] >> n; t1[1] = t0[1] >> n;
		tmp =  t0[1] << (32-n);
		t1[0] |=  tmp;
	}
	else{
		tmp =  t0[1] >> (n-32);
		t1[0] |=  tmp;
	}
	return r;
	
	
	//if (n > 0) for( i=0; i<n; i++ ) z /= 2;
	//else       for( i=0; i<n; i++ ) z *= 2;
	//return z;
#else
	return z >> n;
#endif
	}
*/

#define CHECK_DIFF_NN(  target  , nn_mask ) ( !((nn_mask == ( bitand(nn_mask , target) )) || !( bitand(nn_mask , target) )) )

idx_t ishftc( idx_t i, idx_t s, idx_t L ){
	if(s<0) s += L;
	//return  ((1 << L ) - 1) & ((i << s) | ( i >> (L-s)) );
	return  bitand((bitshift_l(1 , L ) - 1) , bitor(bitshift_l(i , s) , bitshift_r( i , (L-s)) ));
	
}

idx_t power_of_2( idx_t n){
	//return (ghost_gidx_t)(1) << n;
	return bitshift_l( 1 , n );
	}


idx_t Binomial(idx_t N, ghost_gidx_t k){
	
	if((k> N)||(k< 0)) return 0;
	if((k==N)||(k==0)) return 1;
	
	if( 2*k > N ) k=N-k;
	idx_t i;
	idx_t B=1;
	for (i=1;i<=k;i++) B = (B*((N-k)+i))/i;
	
	return B;
	}

idx_t bitcount(idx_t i){
	
	idx_t count = 0;
	//while (count < i){	count++;
	//			i = i & (i-1);  }
	
	while ( i ){	count += bitand(i,1);
			//i >>= 1;
			i = bitshift_r( i,  1 );
			}
	return count;
	}


int SpinChainSZ( ghost_gidx_t row, ghost_lidx_t *nnz, ghost_gidx_t *cols, void *vals){

	static ghost_gidx_t L   = 0;
	static ghost_gidx_t NUp = 0;
	static double Jz  = 1.;
	static double Jxy = 1.;
	static ghost_gidx_t useOBC = 0;
	useOBC &= 1;
	
	
	ghost_gidx_t Ns = Binomial( L, NUp);
	static idx_t * lutlead;
	static idx_t * luttrail;
	static idx_t * cnt;
	static idx_t * lutrev;
	static idx_t * revcnt;
	//static idx_t * imask;
	
	ghost_gidx_t max_row_nnz = 2*L+1;

	if ((row >-1) && (row <Ns)){
		
		double * dvals = vals;
		*nnz = 0;
		
		// binary search
		ghost_gidx_t min_il = 0;
		ghost_gidx_t max_il = power_of_2(L/2);
		ghost_gidx_t il = 0;
		ghost_gidx_t i;
		for(i = 0; i < L+1; i++)
		{
			il = (min_il + max_il)/2;
			if( row >= lutlead[il+1] )
				min_il = il+1;
			else if( max_il != min_il )
				max_il = il;
			else
				break;
		}

		// bad linear search, makes matrix building VERY expensive!
		//  ghost_gidx_t il = 0;
		//  while ( row >= lutlead[il+1] ) il++;
		
		
		idx_t bp;
		//bp = (il << (L/2)) + lutrev[ row - lutlead[il] + revcnt[NUp-bitcount(il)]]; 
		bp = bitshift_l(il , (L/2)) + lutrev[ row - lutlead[il] + revcnt[NUp-bitcount(il)]]; 
		if ( bitcount(bp) != NUp) {  printf("error 1\n"); exit(0);} 
		
		
		ghost_gidx_t hp = 0;
		
		idx_t l;
		
		 if( Jz != 0. ){
			for(l=0;l<L-useOBC;l++){
				//if (bitcount( bp &  ishftc( 1|2 , l , L)) == 1)   hp -= 1;
				if (CHECK_DIFF_NN (bp , ishftc( 1|2 , l , L)) )    hp -= 1;
				else                                              hp += 1;
				
				//int32_t tmp =  ishftc( (bp & imask[l]) , -l , L );   // ( xxXXxxxx & 00110000 ) ->  000000XX
				//switch(tmp){
				//	case 0: hp += 1; break;  // ..00
				//	case 1: hp -= 1; break;  // ..01
				//	case 2: hp -= 1; break;  // ..10
				//	case 3: hp += 1; break;  // ..11
				//	default: printf("error 2\n"); exit(0);
				//}
			}
		 }
		
		if( hp != 0 ){
			cols[ *nnz] = row;
			dvals[*nnz] = hp*Jz;
			*nnz += 1;
		}
		
		if( Jz != 0. ){
			for(l=0;l<L-useOBC;l++){
				
				//int32_t tmp = ishftc( (bp & imask[l]), -l , L);
				//if( (tmp == 1) || (tmp == 2) ){//}
				//if( bitcount( bp &  ishftc( 1|2 , l , L) ) == 1 ){//}
				if (CHECK_DIFF_NN (bp , ishftc( 1|2 , l , L))) {
					
					//int32_t bp2 =  bp ^ imask[l];
					idx_t bp2 =  bitxor(bp , ishftc( 1|2 , l , L));
					
					//ghost_gidx_t lbp = bp2 >> (L/2);
					idx_t lbp = bitshift_r( bp2,  L/2 );
					idx_t tbp = bitand(bp2 , (power_of_2(L/2)-1));
					
					cols[ *nnz] = lutlead[lbp] + luttrail[tbp] -1;
					
					dvals[*nnz] = 2.*Jxy;
					*nnz += 1;
				}
			}
		 }
		
		cols_error_check( Ns , row, *nnz, cols );
		//cols_qsort( cols, vals , sizeof(double) ,  *nnz );
    return 0;
		
	}else if( row == -1 ){
		
		matfuncs_info_t *info = vals;

		info->version   = 1;
		info->base      = 0;
		info->symmetry  = GHOST_SPARSEMAT_SYMM_GENERAL;
		info->datatype  = GHOST_DT_DOUBLE|GHOST_DT_REAL;
		info->nrows     = Ns;
		info->ncols     = Ns;
		info->row_nnz   = max_row_nnz;

		info->hermit    =  1;
		info->eig_info  =  0;
		//info->eig_down  = -8.;
		//info->eig_top   =  8.;

		return 0;
		}

	// =================================================================
	else if ( row == -2 ) {
		L      = cols[0];
		NUp    = cols[1];
		useOBC = cols[2];
		if(L&1) { printf("abort -> need even L!\n"); exit(0); }
		lutlead  = malloc( sizeof(idx_t)* (power_of_2(L/2)+1));
		luttrail = malloc( sizeof(idx_t)*  power_of_2(L/2)   );
		lutrev   = malloc( sizeof(idx_t)* (power_of_2(L/2)+1));
		   cnt   = malloc( sizeof(idx_t)* (L+1) );
		revcnt   = malloc( sizeof(idx_t)* (L+1) );
		//imask    = malloc( sizeof(ghost_gidx_t)*  L    );
		
		*nnz = Binomial( L, NUp);
		
		idx_t i;
		
		for(i=0;i<=L;i++)    cnt[i] = 0;
		for(i=0;i<=L;i++) revcnt[i] = 0;
		
		
		for(i=0;i<power_of_2(L/2);i++){
			idx_t tmp = bitcount(i);
			cnt[tmp]++;
			luttrail[i] = cnt[tmp];
		}
		
		                  revcnt[0] = 0;
		for(i=1;i<=L;i++)  revcnt[i] = revcnt[i-1] + cnt[i-1];
		
		for(i=0;i<power_of_2(L/2);i++){
			idx_t tmp = bitcount(i);
			lutrev[revcnt[tmp]]=i;
			revcnt[tmp]++;
		}
		
		                  revcnt[0] = 0;
		for(i=1;i<=L;i++)  revcnt[i] = revcnt[i-1] + cnt[i-1];
		
		lutlead[0]=0;
		for(i=0;i<power_of_2(L/2);i++){
			idx_t tmp = bitcount(i);
			tmp = NUp-tmp;
			if ( (tmp>=0) && (tmp<=NUp) )   lutlead[i+1] = lutlead[i]+cnt[tmp];
			else                            lutlead[i+1] = lutlead[i];
		}
		
		//for(i=0;i<L-useOBC;i++) imask[i]= ishftc( 1|2 , i, L);  //00000011 00000110, 00001100 .. 11000000, (1000001)
#ifdef DEBUG
        Ns = *nnz;
		printf("N  %ld\n", Ns);
		printf("kl %ld\n", lutlead[power_of_2(L/2)]);
		
		printf("lutlead\n");
		for(i=power_of_2(L/2)-4;i<=power_of_2(L/2);i++)  printf("%ld  ", lutlead[i]);
		
		printf("\n\nluttrail\n");
		for(i=power_of_2(L/2)-4;i< power_of_2(L/2);i++)  printf("%ld  ", luttrail[i]);
		
		printf("\n\nlutrev\n");
		for(i=power_of_2(L/2)-4;i<=power_of_2(L/2);i++)  printf("%ld  ", lutrev[i]);
		
		printf("\n\ncnt\n");
		for(i=0;i<=L;i++)  printf("%ld  ", cnt[i]);
		
		printf("\n\nrevcnt\n");
		for(i=0;i<L;i++)  printf("%ld  ", revcnt[i]);
		
		//printf("\n\nimask\n");
		//for(i=0;i<L;i++)  printf("%d  ", imask[i]);
		
		printf("\n");
#endif
    return 0;
	}
	
	else if( row == -3 ){
		free(lutlead);
		free(luttrail);
		free(lutrev);
		free(cnt);
		free(revcnt);
		//free(imask);
    return 0;
	}

	printf("SpinChainSZ(): error in row %"PRGIDX"\n",row);
	*nnz = -1;
	return 1;              //  error
}



int crsSpinChain( ghost_gidx_t row, ghost_lidx_t *nnz, ghost_gidx_t *cols, void *vals ){

	static ghost_gidx_t L   = 1;
	static double Jx = 1.;
	static double Jy = 1.;
	static double Jz = 1.;
	static double Bx = 1.;
	static double Bz = 1.;
	static int32_t useOBC = 0;
	useOBC &= 1;
	
	ghost_gidx_t  Ns = power_of_2( L );
	
	
	ghost_gidx_t max_row_nnz = 2*L+1;

	if ((row >-1) && (row <Ns)){
		
		double * dvals = vals;
		*nnz = 0;
		ghost_gidx_t il = row;
		
		
		ghost_gidx_t hp_Bz = 0;
		ghost_gidx_t hp_Jz = 0;
		
		int32_t l;
		 if( Bz != 0. ) for(l=0;l<L;l++)if ( (1 << l) &  il ) hp_Bz += 1;
		                                 else                  hp_Bz -= 1;
		
		if( Jz != 0. ){
			for(l=0;l<L-useOBC;l++){
				//if (bitcount( il &  ishftc( 1|2 , l , L)) == 1)   hp_Jz -= 1;
				if ( CHECK_DIFF_NN (il , ishftc( 1|2 , l , L)) ) hp_Jz -= 1;
				else                                              hp_Jz += 1;
			}
		 }
		
		if ( hp_Jz*Jz + hp_Bz*Bz != 0. ){
			cols[ *nnz] = row;
			dvals[*nnz] = hp_Jz*Jz + hp_Bz*Bz;
			*nnz += 1;
		}
		
		if( Bx != 0. ) for(l=0;l<L;l++){
			cols[ *nnz] = il ^ (1 << l);
			dvals[*nnz] = Bx;
			*nnz += 1;
		}
		
		
		if( (Jx != 0.) || (Jy != 0.) ){
			double hp;
			for(l=0;l<L-useOBC;l++){
				//if (bitcount( il &  ishftc( 1|2 , l , L)) == 1)   hp = Jx+Jy;
				if ( CHECK_DIFF_NN (il , ishftc( 1|2 , l , L)) ) hp = Jx+Jy;
				else                                              hp = Jx-Jy;
				
				if ( hp != 0. ){
					cols[ *nnz] = il ^ ishftc( 1|2 , l , L);
					dvals[*nnz] = hp;
					*nnz += 1;
				}
			}
		}
		
		//cols_qsort( cols, vals , sizeof(double) ,  *nnz );
		return 0;
		
	}
	else if( row == -1 ){
		
		matfuncs_info_t *info = vals;

		info->version   = 1;
		info->base      = 0;
		info->symmetry  = GHOST_SPARSEMAT_SYMM_GENERAL;
		info->datatype  = GHOST_DT_DOUBLE|GHOST_DT_REAL;
		info->nrows     = Ns;
		info->ncols     = Ns;
		info->row_nnz   = max_row_nnz;

		info->hermit    =  1;
		info->eig_info  =  0;
		//info->eig_down  = -8.;
		//info->eig_top   =  8.;
		

		return 0;
		}

	// =================================================================
	else if ( row == -2 ) {
		L      = cols[0];
		useOBC = cols[2];
		
		*nnz = power_of_2( L );
		

		return 0;
	}
	
	printf("SpinChainSZ(): error in row %"PRGIDX"\n",row);
	*nnz = -1;
	return 1;              //  error
}

