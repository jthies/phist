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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! compute (shift_j I - A)_i*x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef KACZ_RC_VARIANT
    tmp_r(1:NVEC) = shift_r*x_r(1:NVEC,i)
#else
    tmp_r(1:NVEC) = shift_r*x_r(1:NVEC,i) - &
                    shift_i*x_i(1:NVEC,i)
    tmp_i(1:NVEC) = shift_i*x_r(1:NVEC,i) + &
                    shift_r*x_i(1:NVEC,i)
#endif
#ifdef KACZ_NO_SHIFT
    row_norm=0.0_8
#endif
    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp_r(1:NVEC) = tmp_r(1:NVEC) - val(j)*x_r(1:NVEC,col_idx(j))
#ifdef KACZ_RC_VARIANT
      tmp_i(1:NVEC) = tmp_i(1:NVEC) - val(j)*x_i(1:NVEC,col_idx(j))
#endif
#ifdef KACZ_NO_SHIFT
      row_norm=row_norm+val(j)*val(j)
#endif
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp_r(1:NVEC) = tmp_r(1:NVEC) - val(j)*halo_r(1:NVEC,col_idx(j))
#ifdef KACZ_RC_VARIANT
      tmp_i(1:NVEC) = tmp_i(1:NVEC) - val(j)*halo_i(1:NVEC,col_idx(j))
#endif
#ifdef KACZ_NO_SHIFT
      row_norm=row_norm+val(j)*val(j)
#endif
    end do

#ifndef KACZ_BZERO
    tmp_r(1:NVEC)=tmp_r(1:NVEC)-b(1:NVEC,i)
#endif
    ! Kaczmarz update of X

    ! a) scaling factors
#ifdef KACZ_NO_SHIFT
    tmp_r(1:NVEC)=tmp_r(1:NVEC)*omega/row_norm
#else
    tmp_r(1:NVEC)=tmp_r(1:NVEC)*omega*nrms_ai2i(i)
#endif
#ifdef KACZ_RC_VARIANT
    tmp_i(1:NVEC)=tmp_i(1:NVEC)*omega*nrms_ai2i(i)
#endif
    ! b) projection step
#ifdef KACZ_RC_VARIANT
    x_r(1:NVEC,i)=x_r(1:NVEC,i) - (tmp_r(:)*shift_r+tmp_i(:)*shift_i)
    x_i(1:NVEC,i)=x_i(1:NVEC,i) - (tmp_i(:)*shift_r-tmp_r(:)*shift_i)
#else
    x_r(1:NVEC,i)=x_r(1:NVEC,i) - (tmp_r(:)*shift_r)
#endif
    do j = row_ptr(i), halo_ptr(i)-1, 1
      x_r(1:NVEC,col_idx(j)) = x_r(1:NVEC,col_idx(j)) + &
                               tmp_r(:)*val(j)
#ifdef KACZ_RC_VARIANT
      x_i(1:NVEC,col_idx(j)) = x_i(1:NVEC,col_idx(j)) + &
                               tmp_i(:)*val(j)
#endif
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      halo_r(1:NVEC,col_idx(j)) = halo_r(1:NVEC,col_idx(j)) + &
                               tmp_r(:)*val(j)
#ifdef KACZ_RC_VARIANT
      halo_i(1:NVEC,col_idx(j)) = halo_i(1:NVEC,col_idx(j)) + &
                               tmp_i(:)*val(j)
#endif
    end do
  end do

#ifdef NVEC
#undef NVEC
#endif
