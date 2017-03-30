/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
! this file is included with NVEC #defined to certain common values
! like 1, 2, 4, 8 etc. in kacz_kernels.f90. To get a geral-purpose loop,
! it can be #included without NVEC #defined, in which case a local variable
! NVEC is used. At the end of the file, NVEC is #undefined if it was #defined.
!
! if KACZ_BZERO is #defined, the rhs vector is assumed to be 0 and not accessed.

#ifdef F_DEBUG
#ifdef KACZ_BZERO
write(*,*) 'KACZ-LOOP CLR with b=0, nvec=',NVEC
#else
write(*,*) 'KACZ-LOOP CLR, nvec=',NVEC
#endif  
flush(6)
#endif
  ! shared-memory parallel implementation using coloring info in the map
  do ic = istart, iend,istep
! note: it doesn't make sense to try and use the same scheduling as elsewhere
! (to avoid NUMA problems) because the ordering here is diffrent. For now, we
! will stick with 1 MPI process per socket and think about NUMA later (we could
! 'first PHIST_TOUCH' all vectors using the coloring if it is defined in the map, but
! that would infringe the spMVM performance...)
!$omp parallel do private(tmp_r,tmp_i,i,j,d,row_norm) schedule(static)
  do jc = map%color_offset(ic)+j0,map%color_offset(ic+istep)+j1,istep
    ! two variants: a) seems to be faster on an Emmy Socket. To test the alternative,
    ! uncomment the call permute_local_matrix line in crsmat_module.f90
    ! a) unpermuted matrix, indirect addressing:
    i=map%color_idx(jc)    
    ! b) matrix so that the nodes of a colors form contiguous blocks
    !i=jc
