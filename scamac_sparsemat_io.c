#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>

#include "scamac_include.h"
#include "scamac_include.h"
#include "scamac_sparsemat_io.h"


/* output in Matrix Market format */
int scamac_sparsemat_io_write_mm(const scamac_sparsemat_st *sm, char * fname) {

  if (!sm || !fname) {
    return SCAMAC_EINVAL;
  }

  if (!sm->nr || !sm->nc || !sm->ne) {//empty matrix
    return SCAMAC_EINVAL;
  }

  FILE * f;
  f=fopen(fname,"w"); // suffix ".mm"

  if (sm->valtype == SCAMAC_VAL_REAL) {
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
  } else {
    fprintf(f, "%%%%MatrixMarket matrix coordinate complex general\n");
  }

  fprintf(f, "%"SCAMACPRIDX" %"SCAMACPRIDX" %"SCAMACPRIDX" \n",sm->nr,sm->nc,sm->ne);

  ScamacIdx i,j;
  if (sm->valtype == SCAMAC_VAL_REAL) {
    for (i=0; i<sm->nr; i++) {
      for (j=sm->rptr[i]; j<sm->rptr[i+1]; j++) {
        // MatrixMarket starts with index (1,1)
        fprintf(f, "%"SCAMACPRIDX" %"SCAMACPRIDX" %20.16g\n", i+1, sm->cind[j]+1, sm->val[j]);
      }
    }
  } else {// SCAMAC_VAL_COMPLEX
    for (i=0; i<sm->nr; i++) {
      for (j=sm->rptr[i]; j<sm->rptr[i+1]; j++) {
        // MatrixMarket starts with index (1,1)
        fprintf(f, "%"SCAMACPRIDX" %"SCAMACPRIDX" %20.16g %20.16g\n", i+1, sm->cind[j]+1, sm->val[2*j], sm->val[2*j+1]);
      }
    }
  }

  fclose(f);

  return 0;

}

/* output in Harwell-Boeing format */
int scamac_sparsemat_io_write_hb(const scamac_sparsemat_st *sm, char * fname) {
  if (!sm || !fname) {
    return SCAMAC_EINVAL;
  }

  if (!sm->nr || !sm->nc || !sm->ne) {//empty matrix
    return SCAMAC_EINVAL;
  }

  FILE * f;
  f=fopen(fname,"w");  // suffix ".rb"

  char TITLE[72+1], KEY[8+1];
  ScamacIdx TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD;
  char TOTCRDs[14+1], PTRCRDs[14+1], INDCRDs[14+1], VALCRDs[14+1], RHSCRDs[14+1];
  char MXTYPE[3];
  ScamacIdx NROW, NCOL, NNZERO, NELTVL;
  char NROWs[14+1], NCOLs[14+1], NNZEROs[14+1], NELTVLs[14+1];
  char PTRFMT[16+1], INDFMT[16+1], VALFMT[20+1], RHSFMT[20+1];

  snprintf(TITLE,72,"[SCAMAC] NO TITLE");
  snprintf(KEY,8,"NO KEY");

  if (sm->valtype == SCAMAC_VAL_REAL) {
    MXTYPE[0]='R';
  } else {
    MXTYPE[0]='C';
    printf("%s: Feature not yet implemented\n",__func__);
    exit(EXIT_FAILURE);
  }
  MXTYPE[1]='U'; // we make no use of symmetry
  MXTYPE[2]='A';

  PTRCRD=((sm->nr-1)/8)+1;   // 8 per row
  INDCRD=((sm->ne-1)/8)+1;   // 8 per row
  VALCRD=((sm->ne-1)/3)+1;   // 3 per row
  RHSCRD=0; // no right hand sides
  TOTCRD=PTRCRD+INDCRD+VALCRD+RHSCRD;

  NROW=sm->nr;
  NCOL=sm->nc;
  NNZERO=sm->ne;
  NELTVL=0;  // no elemental matrix


  snprintf(TOTCRDs,14+1,"%"SCAMACPRIDX,TOTCRD);
  snprintf(PTRCRDs,14+1,"%"SCAMACPRIDX,PTRCRD);
  snprintf(INDCRDs,14+1,"%"SCAMACPRIDX,INDCRD);
  snprintf(VALCRDs,14+1,"%"SCAMACPRIDX,VALCRD);
  snprintf(RHSCRDs,14+1,"%"SCAMACPRIDX,RHSCRD);

  snprintf(NROWs,14+1,"%"SCAMACPRIDX,NROW);
  snprintf(NCOLs,14+1,"%"SCAMACPRIDX,NCOL);
  snprintf(NNZEROs,14+1,"%"SCAMACPRIDX,NNZERO);
  snprintf(NELTVLs,14+1,"%"SCAMACPRIDX,NELTVL);

  snprintf(PTRFMT,16+1,"(8I10)"); // sufficient for ne < 10^10
  snprintf(INDFMT,16+1,"(8I10)"); // sufficient for nc,nr < 10^10
  snprintf(VALFMT,20+1,"(3E25.17)"); // sufficient for double precision
  snprintf(RHSFMT,20+1,"(3E25.17)"); // not needed


  // header: 4 lines.
  fprintf(f,"%-72.72s%-8.8s\n",TITLE,KEY);
  fprintf(f,"%-14.14s%-14.14s%-14.14s%-14.14s%-14.14s\n",TOTCRDs,PTRCRDs,INDCRDs,VALCRDs,RHSCRDs);
  fprintf(f,"%-3.3s           %-14.14s%-14.14s%-14.14s%-14.14s\n",MXTYPE,NROWs,NCOLs,NNZEROs,NELTVLs);
  fprintf(f,"%-16.16s%-16.16s%-20.20s%-20.20s\n",PTRFMT,INDFMT,VALFMT,RHSFMT);


  ScamacIdx i,j;
  char OLINE[80+1];
  // data: row pointers
  j=0;
  for (i=0; i<=sm->nr-8; i+=8) { // 8 per row
    // indices start at i
    snprintf(OLINE,80+1,"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX,
             sm->rptr[i]+1, sm->rptr[i+1]+1, sm->rptr[i+2]+1, sm->rptr[i+3]+1, sm->rptr[i+4]+1, sm->rptr[i+5]+1, sm->rptr[i+6]+1, sm->rptr[i+7]+1 );
    fprintf(f,"%s\n",OLINE);
    j+=8; // 8 per row
  }
  // remaining elements
  if (j<=sm->nr) {
    snprintf(OLINE,80+1,"%80.80s"," ");
    for (i=j; i<=sm->nr; i++) {
      snprintf(&OLINE[10*(i-j)],10+1,"%10"SCAMACPRIDX, sm->rptr[i]+1);
    }
    fprintf(f,"%s\n",OLINE);
  }
  // data: column indices
  j=0;
  for (i=0; i<sm->ne-8; i+=8) { // 8 per row
    // indices start at i
    snprintf(OLINE,80+1,"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX"%10"SCAMACPRIDX,
             sm->cind[i]+1, sm->cind[i+1]+1, sm->cind[i+2]+1, sm->cind[i+3]+1, sm->cind[i+4]+1, sm->cind[i+5]+1, sm->cind[i+6]+1, sm->cind[i+7]+1 );
    fprintf(f,"%s\n",OLINE);
    j+=8; // 8 per row
  }
  // remaining elements
  if (j<sm->ne) {
    snprintf(OLINE,80+1,"%80.80s"," ");
    for (i=j; i<sm->ne; i++) {
      snprintf(&OLINE[10*(i-j)],10+1,"%10"SCAMACPRIDX, sm->cind[i]+1);
    }
    fprintf(f,"%s\n",OLINE);
  }
  // data: values
  j=0;
  for (i=0; i<sm->ne-3; i+=3) { // 3 per row
    snprintf(OLINE,80+1,"%25.17e%25.17e%25.17e",
             sm->val[i], sm->val[i+1], sm->val[i+2] );
    fprintf(f,"%s\n",OLINE);
    j+=3; // 3 per row
  }
  // remaining elements
  if (j<sm->ne) {
    snprintf(OLINE,80+1,"%80.80s"," ");
    for (i=j; i<sm->ne; i++) {
      snprintf(&OLINE[25*(i-j)],25+1,"%25.17e", sm->val[i]);
    }
    fprintf(f,"%s\n",OLINE);
  }

  fclose(f);

  return 0;

}

/* output in MATLAB binary format */
int scamac_sparsemat_io_write_matlab(const scamac_sparsemat_st *sm, char * fname) {
  assert(CHAR_BIT == 8);
  if (!sm || !fname) {
    return SCAMAC_EINVAL;
  }

  assert(sizeof *(sm->rptr) == 4);
  assert(sizeof *(sm->cind) == 4);
  assert(sizeof *(sm->val ) == 8);

  if (!sm->nr || !sm->nc || !sm->ne) {//empty matrix
    return SCAMAC_EINVAL;
  }

  FILE * f;
  f=fopen(fname,"wb");

  char buf[256];

  snprintf(buf,128,"%-127.127s","MATLAB 5.0 MAT-file");

  /*
    *( (int16_t *) &buf[124]) =(int16_t) 256;
    *( (int16_t *) &buf[126]) =(int16_t) 19785; // MI
  */

  fwrite(buf, sizeof(char), 124, f);

  int16_t kj[2];
  kj[0]=256;
  kj[1]=19785; // MI
  fwrite(kj, sizeof(int16_t), 2, f);


  int32_t li[4];


  // int ntotbyte = 120;
  int ntotbyte;
  ntotbyte = 16 + 16 + 8 + 3*8; // assuming 4 byte name ['SMAT']
  ntotbyte += 4 * sm->ne + 4 * (sm->nc+1) + 8 * sm->ne;
  if (sm->ne % 2) {// padding for row indices
    ntotbyte +=4;
  }
  if ((sm->nc+1) % 2) {// padding for column indices
    ntotbyte +=4;
  }


  /* sparse matrix (CSC) */
  li[0]=14; // miMATRIX
  li[1]=ntotbyte; // number of bytes
  fwrite(li, sizeof(*li), 2, f);

  /* number of non-zeros */
  li[0]=6; // miUINT32
  li[1]=8; // 8 bytes for field
  li[2]=5 + 256*0; // mxSPARSE_CLASS (real)
  li[3]=sm->ne   ; // ne
  fwrite(li, sizeof(*li), 4, f);

  /* dimension */
  li[0]=5; // miINT32
  li[1]=8; // 8 bytes for field
  li[2]=sm->nr; // dimension
  li[3]=sm->nc; // dimension
  fwrite(li, sizeof(*li), 4, f);

  /* name */
  /*
  li[0]=1 + 256*256*4; // 4 character name
  li[1]=19785;
  fwrite(li, sizeof(*li), 2, f);
  */
  buf[0]=1;
  buf[1]=0;
  buf[2]=4;
  buf[3]=0;
  snprintf(&buf[4],5,"%-4.4s","SMAT");
  fwrite(buf, sizeof(*buf), 8, f);

  /* row indices */
  li[0]=5; // miINT32
  if (sm->ne % 2) {// with padding
    li[1] = 4 * (sm->ne + 1);
    fwrite(li, sizeof(int32_t), 2, f);
    fwrite(sm->cind, sizeof(int32_t), sm->ne-1, f);
    li[0]=sm->cind[sm->ne-1];
    li[1]=0;
    fwrite(li, sizeof(int32_t), 2, f);
  } else {// without padding
    li[1] = 4 * sm->ne;
    fwrite(li, sizeof(int32_t), 2, f);
    fwrite(sm->cind, sizeof(int32_t), sm->ne, f);
  }

  /* column indices */
  li[0]=5; // miINT32
  if ((sm->nc+1) % 2) {// with padding
    li[1] = 4 * (sm->nc + 2);
    fwrite(li, sizeof(int32_t), 2, f);
    fwrite(sm->rptr, sizeof(int32_t), sm->nc, f);
    li[0]=sm->rptr[sm->nc];
    li[1]=0;
    fwrite(li, sizeof(int32_t), 2, f);
  } else {// without padding
    li[1] = 4 * (sm->nc + 1);
    fwrite(li, sizeof(int32_t), 2, f);
    fwrite(sm->rptr, sizeof(int32_t), sm->nc+1, f);
  }

  /* values */
  li[0]=9; // miDOUBLE
  li[1]=8 * sm->ne; // number of bytes
  fwrite(li, sizeof(*li), 2, f);
  fwrite(sm->val, sizeof(double), sm->ne, f);

  fclose(f);

  return 0;

}

/* output in GHOST binary format */
int scamac_sparsemat_io_write_ghost(const scamac_sparsemat_st *sm, char * fname) {
  if (!sm || !fname) {
    return SCAMAC_EINVAL;
  }

  if (!sm->nr || !sm->nc || !sm->ne) {//empty matrix
    return SCAMAC_EINVAL;
  }

  if ( sizeof *(sm->rptr) != 8 || sizeof *(sm->cind) != 8 || sizeof *(sm->val) !=8) {
    printf("%s: Wrong size of sparse matrix. Abort.\n",__func__);
    exit(EXIT_FAILURE);
  }

  FILE * f;
  f=fopen(fname,"wb");  // suffix ".ghm"

  // write header

  int16_t shint[5];
  // Endian indicator (here: Little Endian assumed)
  shint[0]=0;
  // version data format
  shint[1]=1;
  // C: Start counting with 0
  shint[2]=0;
  // general matrix (no symmetry)
  shint[3]=1;
  // data type (here: real double)
  shint[4]=6;
  fwrite(shint, sizeof (int16_t), 5, f);

  int64_t llint[3];
  llint[0]=sm->nr;
  llint[1]=sm->nc;
  llint[2]=sm->ne;
  fwrite(llint, sizeof (int64_t), 3, f);

  // write data

  //row pointers
  fwrite(sm->rptr, sizeof *(sm->rptr), sm->nr+1, f);
  //column indices
  fwrite(sm->cind, sizeof *(sm->cind), sm->ne, f);
  //values
  fwrite(sm->val, sizeof *(sm->val), sm->ne, f);

  fclose(f);

  return 0;

}

