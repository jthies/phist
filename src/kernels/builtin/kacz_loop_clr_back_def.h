! this file is included with NVEC #defined to certain common values
! like 1, 2, 4, 8 etc. in kacz_kernels.f90. To get a geral-purpose loop,
! it can be #included without NVEC #defined, in which case a local variable
! NVEC is used. At the end of the file, NVEC is #undefined if it was #defined.
!
! if KACZ_BZERO is #defined, the rhs vector is assumed to be 0 and not accessed.

#ifdef F_DEBUG
#ifdef KACZ_BZERO
write(*,*) 'KACZ-LOOP CLR back with b=0, nvec=',NVEC
#else
write(*,*) 'KACZ-LOOP CLR back, nvec=',NVEC
#endif  
flush(6)
#endif
  ! shared-memory parallel implementation using coloring info in the map
  do ic = map%ncolors, 1,-1
! note: it doesn't make sense to try and use the same scheduling as elsewhere
! (to avoid NUMA problems) because the ordering here is diffrent. For now, we
! will stick with 1 MPI process per socket and think about NUMA later (we could
! 'first PHIST_TOUCH' all vectors using the coloring if it is defined in the map, but
! that would infringe the spMVM performance...)
!$omp parallel do private(tmp_r,tmp_i,i,j,d,row_norm) schedule(static)
  do jc = map%color_offset(ic),map%color_offset(ic+1)-1,1
    i=map%color_idx(jc)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! compute (shift_j I - A)_i*x
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    d=0.0_8
    row_norm(1:NVEC)=0.0_8 ! compute row 2-norm on-the-fly
#ifdef KACZ_NO_SHIFT
    tmp_r(1:NVEC) = 0.0_8
# ifdef KACZ_RC_VARIANT
    tmp_i(1:NVEC) = 0.0_8
# endif
#else
# ifndef KACZ_RC_VARIANT
    tmp_r(1:NVEC) = shift_r(1:NVEC)*x_r(1:NVEC,i)
# else
    tmp_r(1:NVEC) = shift_r(1:NVEC)*x_r(1:NVEC,i) - &
                    shift_i(1:NVEC)*x_i(1:NVEC,i)
    tmp_i(1:NVEC) = shift_i(1:NVEC)*x_r(1:NVEC,i) + &
                    shift_r(1:NVEC)*x_i(1:NVEC,i)
# endif
#endif

    do j = row_ptr(i), halo_ptr(i)-1, 1
      tmp_r(1:NVEC) = tmp_r(1:NVEC) - val(j)*x_r(1:NVEC,col_idx(j))
#ifdef KACZ_RC_VARIANT
      tmp_i(1:NVEC) = tmp_i(1:NVEC) - val(j)*x_i(1:NVEC,col_idx(j))
#endif
      row_norm(1)=row_norm(1)+val(j)*val(j)
#ifndef KACZ_NO_SHIFT
      if (col_idx(j)==j) then
        d=val(j)
      end if
#endif
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      tmp_r(1:NVEC) = tmp_r(1:NVEC) - val(j)*halo_r(1:NVEC,col_idx(j))
#ifdef KACZ_RC_VARIANT
      tmp_i(1:NVEC) = tmp_i(1:NVEC) - val(j)*halo_i(1:NVEC,col_idx(j))
#endif
      row_norm(1)=row_norm(1)+val(j)*val(j)
    end do

#ifndef KACZ_BZERO
    tmp_r(1:NVEC)=tmp_r(1:NVEC)-b(1:NVEC,i)
#endif

#ifdef KACZ_NO_SHIFT
  row_norm(1:NVEC) = row_norm(1)
#else
  ! correct for real or complex shift
  row_norm(1:NVEC) = row_norm(1) + 2*d*shift_r(1:NVEC)+shift_r(1:NVEC)*shift_r(1:NVEC)
# ifdef KACZ_RC_VARIANT
  row_norm(1:NVEC) = row_norm(1:NVEC)+shift_i(1:NVEC)*shift_i(1:NVEC)
# endif
#endif

    ! Kaczmarz update of X

    ! a) scaling factors
    tmp_r(1:NVEC)=tmp_r(1:NVEC)*omega(1:NVEC)/row_norm(1:NVEC)
#ifdef KACZ_RC_VARIANT
    tmp_i(1:NVEC)=tmp_i(1:NVEC)*omega(1:NVEC)/row_norm(1:NVEC)
#endif
    ! b) projection step
#ifndef KACZ_NO_SHIFT
# ifdef KACZ_RC_VARIANT
    x_r(1:NVEC,i)=x_r(1:NVEC,i) - (tmp_r(1:NVEC)*shift_r+tmp_i(1:NVEC)*shift_i)
    x_i(1:NVEC,i)=x_i(1:NVEC,i) - (tmp_i(1:NVEC)*shift_r-tmp_r(1:NVEC)*shift_i)
# else
    x_r(1:NVEC,i)=x_r(1:NVEC,i) - (tmp_r(1:NVEC)*shift_r)
# endif
#endif
    do j = row_ptr(i), halo_ptr(i)-1, 1
      x_r(1:NVEC,col_idx(j)) = x_r(1:NVEC,col_idx(j)) + &
                               tmp_r(1:NVEC)*val(j)
#ifdef KACZ_RC_VARIANT
      x_i(1:NVEC,col_idx(j)) = x_i(1:NVEC,col_idx(j)) + &
                               tmp_i(1:NVEC)*val(j)
#endif
    end do
    do j = halo_ptr(i), row_ptr(i+1)-1, 1
      halo_r(1:NVEC,col_idx(j)) = halo_r(1:NVEC,col_idx(j)) + &
                               tmp_r(1:NVEC)*val(j)
#ifdef KACZ_RC_VARIANT
      halo_i(1:NVEC,col_idx(j)) = halo_i(1:NVEC,col_idx(j)) + &
                               tmp_i(1:NVEC)*val(j)
#endif
    end do
  end do
  end do

#ifdef NVEC
#undef NVEC
#endif
