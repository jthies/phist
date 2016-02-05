#include "matfuncs.h"
#include "part_tools.h"

//#define WRITE_MATRIX
#define REPARTITION

#ifndef REPARTITION
# define PERM( _row ) _row
# define IPERM( _row ) _row
#else
# define PERM( _row ) perm2d(_row, +1)
# define IPERM( _row ) perm2d(_row, -1)
#endif

#ifndef graphene_a_unit
#define graphene_a_unit 0.142 //  carbon-carbon bond length in nanometers
#endif

#ifndef SQRT3
#define SQRT3 1.732050807568877
#endif

#ifndef D_ZZ
#define D_ZZ  (0.5*graphene_a_unit*SQRT3)
#endif
#ifndef D_AC
#define D_AC  (1.5*graphene_a_unit)
#endif

int crsGraphene( ghost_gidx row, ghost_lidx *nnz, ghost_gidx *cols, void *vals, void* data){

	static ghost_gidx L = 4 ;
	static ghost_gidx W = 4 ;
	ghost_gidx N = L*W;

	static int32_t zigzag_first         = 1;
	static int32_t long_range_hopping   = 1;
	static int32_t PBC_W                = 1;
	static int32_t PBC_L                = 1;
	static int32_t basis_place          = 0;
	static int32_t use_symmflag         = 0;
	static int32_t J_mode           = 0;
	
	static	double gamma = 0.;  // disorder
	
  //static	double Dot_R = 0.;
  //static	double Dot_V = 0.;
  //static int32_t use_dots         = 0;
	
	static double t1 = -1.;
	static double t2 = -0.1/2.78;;
	static double t3 = -0.095/2.78;

	ghost_gidx           max_row_nnz  = 4 ;
	if(long_range_hopping) max_row_nnz += 9 ;


	if ((row >-1 ) && (row <N)){       //  defined output -- write entries of #row in *cols and *vals
	                                   //                    return number of entries
		double * dvals = vals;
		ghost_gidx i = 0 ;
		row = IPERM(row);


		ghost_gidx w = row%W;
		ghost_gidx l = row/W;

		int32_t   uplink = (int32_t)( (l+w+basis_place  )& 1  );
		
		
		double V = gamma*( ((double)(rand()))/((double)(RAND_MAX)) -0.5);
		//if ( use_dots ){
		//double Dot_D = 4*Dot_R;
		//double Dot_cx= 2*Dot_R;
		//double Dot_cy= 2*Dot_R;
		//double P_x,P_y;
		//if( zigzag_first ) graphene_mesh_coo( w , l, basis_place, &P_x, &P_y);
		//else               graphene_mesh_coo( l , w, basis_place, &P_x, &P_y);
		
		//double Dot_x = fmod( P_x , Dot_D );
		//double Dot_y = fmod( P_y , Dot_D );
		
		//if( (((Dot_cx-Dot_x)*(Dot_cx-Dot_x) + (Dot_cy-Dot_y)*(Dot_cy-Dot_y)) < Dot_R*Dot_R )  )  V += Dot_V;
		//}
		
		double fac_Jzz[5] = {1.,1.,1.,1.,1.};
		double fac_Jac[9]= {1.,1.,1.,1.,1.,1.,1.,1.,1.};
		if(J_mode & 1) {fac_Jzz[0] = -SQRT3; fac_Jzz[1] = -SQRT3*0.5; fac_Jzz[2] = 0.; fac_Jzz[3] = SQRT3*0.5; fac_Jzz[4] = SQRT3;}
		if(J_mode & 2) { fac_Jac[0] = -2.0;
		                 fac_Jac[1] = -1.5;
		                 fac_Jac[2] = -1.0;
		                 fac_Jac[3] = -0.5;
		                 fac_Jac[4] =  0.0;
		                 fac_Jac[5] =  0.5;
		                 fac_Jac[6] =  1.0;
		                 fac_Jac[7] =  1.5;
		                 fac_Jac[8] =  2.0; }
		

		if( zigzag_first ){
			if(!use_symmflag){
			    if( PBC_L && (l==L-1 ))      {
			      if( long_range_hopping ){
			        if(PBC_W && (w >W-2 ) )            {   cols[i] = row + W-N + 1 -W;   dvals[i] = t2 *fac_Jac[7]*fac_Jzz[3]; i++; }
			        if(PBC_W && (w >W-3 ) && uplink )  {   cols[i] = row + W-N + 2 -W;   dvals[i] = t3 *fac_Jac[6]*fac_Jzz[4]; i++; }
			        if((w>1) && uplink )               {   cols[i] = row + W-N - 2 ;     dvals[i] = t3 *fac_Jac[6]*fac_Jzz[0]; i++; }
			        if( w>0  )                         {   cols[i] = row + W-N - 1 ;     dvals[i] = t2 *fac_Jac[7]*fac_Jzz[1]; i++; }
			      }
			      if( uplink )                         {   cols[i] = row + W-N;          dvals[i] = t1 *fac_Jac[6]*fac_Jzz[2]; i++; }
			      if(!uplink && long_range_hopping )   {   cols[i] = row + W-N;          dvals[i] = t3 *fac_Jac[8]*fac_Jzz[2]; i++; }
			      if( long_range_hopping ){
			        if( w <W-1  )                      {   cols[i] = row + W-N + 1 ;     dvals[i] = t2 *fac_Jac[7]*fac_Jzz[3]; i++; }
			        if( (w <W-2) && uplink )           {   cols[i] = row + W-N + 2 ;     dvals[i] = t3 *fac_Jac[6]*fac_Jzz[4]; i++; }
			        if(PBC_W &&  (w<2) && uplink  )    {   cols[i] = row + W-N - 2 +W;   dvals[i] = t3 *fac_Jac[6]*fac_Jzz[0]; i++; }
			        if(PBC_W &&  (w<1) )               {   cols[i] = row + W-N - 1 +W;   dvals[i] = t2 *fac_Jac[7]*fac_Jzz[2]; i++; }
			      }
			      }
			    //---------------------------------------------------------------------------------------
			    if ( l>0  ) {
			      if( long_range_hopping ){
			        if(PBC_W && (w >W-2 ) )            {   cols[i] = row - W + 1 -W;     dvals[i] = t2 *fac_Jac[1]*fac_Jzz[3]; i++; }
			        if(PBC_W && (w >W-3 ) && !uplink ) {   cols[i] = row - W + 2 -W;     dvals[i] = t3 *fac_Jac[2]*fac_Jzz[4]; i++; }
			        if((w>1) && !uplink  )             {   cols[i] = row - W - 2 ;       dvals[i] = t3 *fac_Jac[2]*fac_Jzz[0]; i++; }
			        if( w>0  )                         {   cols[i] = row - W - 1 ;       dvals[i] = t2 *fac_Jac[1]*fac_Jzz[1]; i++; }
			      }
			     if(!uplink )                          {   cols[i] = row - W;            dvals[i] = t1 *fac_Jac[2]*fac_Jzz[2]; i++; }
			     if( uplink && long_range_hopping )    {   cols[i] = row - W;            dvals[i] = t3 *fac_Jac[0]*fac_Jzz[2]; i++; }
			      if( long_range_hopping ){
			        if( w <W-1  )                      {   cols[i] = row - W + 1 ;       dvals[i] = t2 *fac_Jac[1]*fac_Jzz[3]; i++; }
			        if((w <W-2) && !uplink )           {   cols[i] = row - W + 2 ;       dvals[i] = t3 *fac_Jac[2]*fac_Jzz[4]; i++; }
			        if(PBC_W &&  (w<2) && !uplink )    {   cols[i] = row - W - 2 +W;     dvals[i] = t3 *fac_Jac[2]*fac_Jzz[0]; i++; }
			        if(PBC_W && (w<1) )                {   cols[i] = row - W - 1 +W;     dvals[i] = t2 *fac_Jac[1]*fac_Jzz[1]; i++; }
			      }
				 }
			    //--------------------------------------------------------------------------------------
			    if(PBC_W){
			      if( w >W-2L )                        {   cols[i] = row + 1 -W;         dvals[i] = t1 *fac_Jac[5-2*uplink]*fac_Jzz[3]; i++; }
			      if((w >W-3L) && long_range_hopping)  {   cols[i] = row + 2 -W;         dvals[i] = t2 *fac_Jac[4]*fac_Jzz[4]; i++; }
			      }
			    if((w>1 ) && long_range_hopping )      {   cols[i] = row - 2 ;           dvals[i] = t2 *fac_Jac[4]*fac_Jzz[0]; i++; }
			    if( w>0  )                             {   cols[i] = row - 1 ;           dvals[i] = t1 *fac_Jac[5-2*uplink]*fac_Jzz[1]; i++; }
			 }
			    //======================================================================================
			    // todo onsite 
			                                               if (long_range_hopping) V += 3*t3;
			                                               cols[i] = row ;               dvals[i] = V *fac_Jac[4]*fac_Jzz[2]; i++; 
			    //======================================================================================
			    if( w <W-1  )                          {   cols[i] = row + 1 ;           dvals[i] = t1 *fac_Jac[5-2*uplink]*fac_Jzz[3]; i++; }
			    if((w <W-2 ) && long_range_hopping )   {   cols[i] = row + 2 ;           dvals[i] = t2 *fac_Jac[4]*fac_Jzz[4]; i++; }
			    if(PBC_W){
			       if((w <2 ) && long_range_hopping  ) {   cols[i] = row - 2 +W;         dvals[i] = t2 *fac_Jac[4]*fac_Jzz[0]; i++; }
			       if( w==0  )                         {   cols[i] = row - 1 +W;         dvals[i] = t1 *fac_Jac[5-2*uplink]*fac_Jzz[1]; i++; }
			      }
			    //--------------------------------------------------------------------------------------
			    if( l<L-1 ){
			      if( long_range_hopping ){
			        if(PBC_W && (w >W-2 ) )            {   cols[i] = row + W + 1 -W;     dvals[i] = t2 *fac_Jac[7]*fac_Jzz[3]; i++; }
			        if(PBC_W && (w >W-3 ) && uplink )  {   cols[i] = row + W + 2 -W;     dvals[i] = t3 *fac_Jac[6]*fac_Jzz[4]; i++; }
			        if((w>1) && uplink )               {   cols[i] = row + W - 2 ;       dvals[i] = t3 *fac_Jac[6]*fac_Jzz[0]; i++; }
			        if( w>0  )                         {   cols[i] = row + W - 1 ;       dvals[i] = t2 *fac_Jac[7]*fac_Jzz[1]; i++; }
			      }
			      if(uplink )                          {   cols[i] = row + W;            dvals[i] = t1 *fac_Jac[6]*fac_Jzz[2]; i++; }
			      if(!uplink && long_range_hopping )   {   cols[i] = row + W;            dvals[i] = t3 *fac_Jac[8]*fac_Jzz[2]; i++; }
			      if( long_range_hopping ){
			        if(  w <W-1  )                     {   cols[i] = row + W + 1 ;       dvals[i] = t2 *fac_Jac[7]*fac_Jzz[3]; i++; }
			        if( (w <W-2) && uplink )           {   cols[i] = row + W + 2 ;       dvals[i] = t3 *fac_Jac[6]*fac_Jzz[4]; i++; }
			        if(PBC_W &&  (w<2) && uplink  )    {   cols[i] = row + W - 2 +W;     dvals[i] = t3 *fac_Jac[6]*fac_Jzz[0]; i++; }
			        if(PBC_W &&  (w<1) )               {   cols[i] = row + W - 1 +W;     dvals[i] = t2 *fac_Jac[7]*fac_Jzz[1]; i++; }
			      }
			     }
			    //---------------------------------------------------------------------------------------
			    if(PBC_L && (l==0 )){
			      if( long_range_hopping ){
			        if(PBC_W && (w >W-2 ) )            {   cols[i] = row - W+N + 1 -W;   dvals[i] = t2 *fac_Jac[1]*fac_Jzz[3]; i++; }
			        if(PBC_W && (w >W-3 ) && !uplink ) {   cols[i] = row - W+N + 2 -W;   dvals[i] = t3 *fac_Jac[2]*fac_Jzz[4]; i++; }
			        if((w>1) && !uplink  )             {   cols[i] = row - W+N - 2 ;     dvals[i] = t3 *fac_Jac[2]*fac_Jzz[0]; i++; }
			        if( w>0  )                         {   cols[i] = row - W+N - 1 ;     dvals[i] = t2 *fac_Jac[1]*fac_Jzz[1]; i++; }
			      }
			      if(!uplink )                         {   cols[i] = row - W+N;          dvals[i] = t1 *fac_Jac[2]*fac_Jzz[2]; i++; }
			      if( uplink && long_range_hopping )   {   cols[i] = row - W+N;          dvals[i] = t3 *fac_Jac[0]*fac_Jzz[2]; i++; }
			      if( long_range_hopping ){
			        if( w <W-1  )                      {   cols[i] = row - W+N + 1 ;     dvals[i] = t2 *fac_Jac[1]*fac_Jzz[3]; i++; }
			        if((w <W-2) && !uplink )           {   cols[i] = row - W+N + 2 ;     dvals[i] = t3 *fac_Jac[2]*fac_Jzz[4]; i++; }
			        if(PBC_W &&  (w<2) && !uplink )    {   cols[i] = row - W+N - 2 +W;   dvals[i] = t3 *fac_Jac[2]*fac_Jzz[0]; i++; }
			        if(PBC_W &&  (w<1) )               {   cols[i] = row - W+N - 1 +W;   dvals[i] = t2 *fac_Jac[1]*fac_Jzz[1]; i++; }
			      }
		            }

			}
		//else{
			// todo  !zigzag_first    }

		ghost_gidx j;
		if( i>max_row_nnz )
      printf("corrupt max_row_nnz for row %" PRGIDX "   nnz = %" PRGIDX " \n",row, i);
		for( j=0;j<i;j++ )
      if( (cols[j]<0) || (cols[j]>=N) )
        printf("corrupt col_idx in row %" PRGIDX "    (w,l) =(%" PRGIDX ",%" PRGIDX ") : col[%" PRGIDX "] = %" PRGIDX "\n", row, w,l, j, cols[j]);

		*nnz = (ghost_lidx)i;
#ifdef WRITE_MATRIX
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	char fname[200];
	sprintf(fname,"mat%d.txt",rank);
	FILE *matfile=fopen(fname,"a");
	PHIST_OUT(PHIST_DEBUG,"open file '%s'\n",fname);
	if (!matfile) PHIST_OUT(PHIST_ERROR,"could not open file\n");
	fprintf(matfile,"%ld:",row);
	for (i=0;i<*nnz;i++) fprintf(matfile," %ld (%ld)",cols[i], PERM(cols[i]));
	fprintf(matfile,"\n");
	fclose(matfile);
#endif
	for (i=0;i<*nnz;i++) cols[i]=PERM(cols[i]);
		return 0;
	}

	// =================================================================
	else if ( row == -1) {
		matfuncs_info_t *info = vals;

		info->version   = 1;
		info->base      = 0;
		info->symmetry  = GHOST_SPARSEMAT_SYMM_GENERAL;
		info->datatype  = GHOST_DT_DOUBLE|GHOST_DT_REAL;
		info->nrows     = N;
		info->ncols     = N;
		info->row_nnz   = max_row_nnz;

		info->hermit    =  1;
		info->eig_info  =  1;
		info->eig_down  = -4.;
		info->eig_top   =  4.;

		return 0;
	}else if ( row == -2) {
		W = nnz[0];
		L = nnz[1];
		perm2d(W,L);

		//printf("W,L = %d, %d\n",W,L);
		//printf("nnz = %d\n",max_row_nnz);
		return 0;
	}else if ( row == -3) {
		double * dvals = vals;
		W     = (ghost_gidx)(dvals[0]/D_ZZ);
		L     = (ghost_gidx)(dvals[1]/D_AC);
		L += L&1;
		W += W&1;
		
    //Dot_R = dvals[2];
    //Dot_V = dvals[3];
    //use_dots = 1;
		
		cols[0] = W;
		cols[1] = L;
		
		*nnz = W*L;
		

		//printf("W,L = %d, %d\n",W,L);
		//printf("nnz = %d\n",max_row_nnz);
		return 0;
	}else if ( row == -4) {
		
		PBC_W = (int32_t)(1 & nnz[0]);
		PBC_L = (int32_t)(1 & nnz[1]);
		
		return 0;
	}else if ( row == -5) {
		
		gamma = *((double *)(vals));
		
		return 0;
	}
	
	

	printf("################### error \n");
	*nnz = -1; //  error
	return 1;
}


