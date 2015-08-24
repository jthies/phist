!> \file random.f90
!! 'parallel' random number generator
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!! code for KISS prng copied from http://web.mst.edu/~vojtat/class_5403/kiss05/rkiss05.f90
!!

#include "phist_config.h"
module random_module
  implicit none
  private

  public :: phist_random_seed, phist_Drandom_number, drandom_1, drandom_general

  integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
  integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

  integer(kind=i4b) :: state_x, state_y, state_z, state_w

contains

!! copied from http://web.mst.edu/~vojtat/class_5403/kiss05/rkiss05.f90
!! (with some modifications)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis  
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123  
! 
! 
! A call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05 call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
! 
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
! 
! 
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit 
!        v093     Aug 13, 2012    changed inter representation test to avoid data statements
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


pure subroutine rkiss05(x,y,z,w, val)
implicit none

real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

real(r8b), intent(out)       :: val
integer(i4b),intent(inout)   :: x,y,z,w              ! working variables for the four generators
integer(i4b)          :: kiss

x = 69069 * x + 1327217885
y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
z = 18000 * iand (z, 65535) + ishft (z, - 16)
w = 30903 * iand (w, 65535) + ishft (w, - 16)
kiss = ishft(x + y + ishft (z, 16) + w , -1)
val=2.0_r8b*(kiss*am-0.5_r8b)
end subroutine rkiss05


! TODO: this can be done in O(log n) time by some more sophisticated methods!
pure subroutine kissskip(n,x,y,z,w)
implicit none

integer(kind=8), intent(in)   :: n
integer(i4b), intent(inout)   :: x,y,z,w              ! working variables for the four generators
integer(kind=8) :: i

do  i = 1, n, 1
  x = 69069 * x + 1327217885
  y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
  z = 18000 * iand (z, 65535) + ishft (z, - 16)
  w = 30903 * iand (w, 65535) + ishft (w, - 16)
end do
end subroutine kissskip


! TODO: initialization with 4 integers!
subroutine kissinit(iinit,x,y,z,w)
implicit none

integer(i4b) idum,ia,im,iq,ir,iinit
integer(i4b) k,x,y,z,w,c1
parameter (ia=16807,im=2147483647,iq=127773,ir=2836)

!!! Test integer representation !!!
c1=-8
c1=ishftc(c1,-3)
!     print *,c1
if (c1.ne.536870911) then
   print *,'Nonstandard integer representation. Stoped.'
   stop
endif

idum=iinit
idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
if (idum.eq.0) idum=1
if (idum.ge.IM) idum=IM-1

k=(idum)/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum = idum + IM
if (idum.lt.1) then
   x=idum+1 
else 
   x=idum
endif
k=(idum)/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum = idum + IM
if (idum.lt.1) then 
   y=idum+1 
else 
   y=idum
endif
k=(idum)/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum = idum + IM
if (idum.lt.1) then
   z=idum+1 
else 
   z=idum
endif
k=(idum)/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum = idum + IM
if (idum.lt.1) then
   w=idum+1 
else 
   w=idum
endif

call kissskip(1_8,x,y,z,w)

return
end subroutine kissinit

! seed random number generator
subroutine phist_random_seed(seed) bind(C)
  use iso_c_binding, only: C_INT
  integer(C_INT), value :: seed
  integer(kind=i4b) :: x, y, z, w

  call kissinit(seed, x, y, z, w)

!$omp critical (random_module_state)
  state_x = x
  state_y = y
  state_z = z
  state_w = w
!$omp end critical (random_module_state)

end subroutine phist_random_seed


subroutine phist_Drandom_number(n, r) bind(C, name='phist_Drandom_number')
  use iso_c_binding, only: C_INT, C_DOUBLE
  integer(kind=C_INT), value :: n
  real(kind=C_DOUBLE) :: r(n)
  integer :: i


!$omp critical (random_module_state)
  do i = 1, n, 1
    call rkiss05(state_x,state_y,state_z,state_w,r(i))
  end do
!$omp end critical (random_module_state)

end subroutine phist_Drandom_number


subroutine drandom_1(nrows, y, pre_skip, post_skip) bind(C,name='drandom_1_')
  use iso_c_binding, only: C_INT, C_INT64_T, C_DOUBLE
  implicit none
  integer(kind=C_INT), intent(in) :: nrows
  integer(kind=C_INT64_T), intent(in) :: pre_skip, post_skip
  real(kind=C_DOUBLE), intent(out) :: y(*)
!dir$ assume_aligned y:64

  call drandom_general(1, nrows, y, 1, pre_skip, post_skip)

end subroutine drandom_1


subroutine drandom_general(nvec, nrows, y, ldy, pre_skip, post_skip) bind(C,name='drandom_general_')
  use iso_c_binding, only: C_INT, C_INT64_T, C_DOUBLE
#ifdef PHIST_HAVE_OPENMP
  use omp_lib
#endif
  implicit none
  integer(kind=C_INT), intent(in) :: nvec, nrows, ldy
  integer(kind=C_INT64_T), intent(in) :: pre_skip, post_skip
  real(kind=C_DOUBLE), intent(out) :: y(ldy,*)
  integer :: i, j, it, nt
  integer(kind=i4b) :: sx, sy, sz, sw
  integer :: rowsPerThread, rowsModThread
  integer, allocatable :: threadOffset(:)
!dir$ assume_aligned y:8

#ifdef PHIST_HAVE_OPENMP
  nt = omp_get_max_threads()
#else
  nt = 1
#endif
  allocate(threadOffset(0:nt))
  rowsPerThread = nrows/nt
  rowsModThread = mod(nrows, nt)
  threadOffset(0) = 0
  do i = 1, rowsModThread, 1
    threadOffset(i) = threadOffset(i-1) + rowsPerThread+1
  end do
  do i = rowsModThread+1, nt, 1
    threadOffset(i) = threadOffset(i-1) + rowsPerThread
  end do



!$omp critical (random_module_state)
  sx = state_x
  sy = state_y
  sz = state_z
  sw = state_w
! handcoded openmp parallel do schedule(static)
!$omp parallel private(i,j,it) firstprivate(sx,sy,sz,sw)
#ifdef PHIST_HAVE_OPENMP
  it = omp_get_thread_num()
#else
  it = 0
#endif

  ! jump to own offset
  call kissskip(pre_skip+threadOffset(it)*nvec,sx,sy,sz,sw)

  do i = threadOffset(it)+1, threadOffset(it+1), 1
    do j = 1, nvec
      call rkiss05(sx,sy,sz,sw,y(j,i))
    end do
  end do

  ! update state on last thread
  if( it .eq. nt-1 ) then
    state_x = sx
    state_y = sy
    state_z = sz
    state_w = sw
  end if
!$omp end parallel

  call kissskip(post_skip,state_x,state_y,state_z,state_w)

!$omp end critical (random_module_state)

end subroutine drandom_general

end module random_module

