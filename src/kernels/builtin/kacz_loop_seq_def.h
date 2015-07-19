! this file is included with NVEC #defined to certain common values
! like 1, 2, 4, 8 etc. in kacz_kernels.f90. To get a geral-purpose loop,
! it can be #included without NVEC #defined, in which case a local variable
! NVEC is used. At the end of the file, NVEC is #undefined if it was #defined.
!
! if KACZ_BZERO is #defined, the rhs vector is assumed to be 0 and not accessed.

#ifdef F_DEBUG
#ifdef KACZ_BZERO
write(*,*) 'KACZ-LOOP SEQ with b=0, nvec=',NVEC
#else
write(*,*) 'KACZ-LOOP SEQ, nvec=',NVEC
#endif  
flush(6)
#endif
  
  ! sequential implementation ignoring coloring information
  ! (lexicographic or given ordering)
  do i = istart, iend,istep
