#include "phist_config.h"
#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#endif

#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
#define G_GIDX_T C_INT64_T
#else
#ifdef GHOST_HAVE_LONGIDX_LOCAL
#define G_LIDX_T C_INT64_T
#else
#define G_LIDX_T C_INT32_T
#endif
#ifdef GHOST_HAVE_LONGIDX_GLOBAL
#define G_GIDX_T C_INT64_T
#else
#define G_GIDX_T C_INT32_T
#endif
#endif

! Generalization of MATPDE to 3D
! (from the NEP matrix collection in the matrix market MATPDE matrix generator)
! We allow generating various test problems by calling some initialization
! routines and provide simple matfuncs in matfuncs.h to wrap this functionality.
!
! The currently implemented test cases are:

!=======================================================================================
! (A) from Gordon & Gordon CARP-CG Paper (Parallel Computing 2009)                      
!                                                                                       
!0. Du                                                                              = F 
!1. Du +          1000u_x                                                           = F 
!2. Du + 1000exp(xyz)(u_x   +               u_y -      u_z)                         = F 
!3. Du +          100xu_x   -              yu_y +     zu_z  + 100(x + y + z)/(xyz)u = F 
!4. Du -      10^5x^2(u_x   +               u_y +      u_z)                         = F 
!5. Du -   1000(1+x^2)u_x   +           100(u_y +      u_z)                         = F 
!6. Du -   1000[(1-2x)u_x   +         (1-2y)u_y+ (1-2z)u_z]                         = F 
!7. Du -       1000x^2u_x                                   +                 1000u = F 
!8. Du - d(  10exp(xy)u)/dx - d(  10exp(-xy)u)/dy                                   = F 
!9. Du - d(1000exp(xy)u)/dx - d(1000exp(-xy)u)/dy                                   = F 
!                                                                                       
!analytic solutions                                                                     
!                                                                                       
!Problem 1: u(x,y,z) = xyz(1-x)(1-y)(1-z)                                               
!Problem 2: u(x,y,z) = x + y + z                                                        
!   "  3–7: u(x,y,z) = exp(xyz)sin(pi x)sin(pi y)sin(pi z)                              
!                                                                                       
!                                                                                       
!boundary conditions: Dirichlet                                                         
!                                                                                       
!Problem 8,9: set u=1, b=A*u and use u=0 as BC                                          
!                                                                                       
!=======================================================================================
! (B) problems with variable coefficentts in the diffusion term                         
! 0. -(exp(-xyz) Ux)x - (exp(xyz) Uy)y  -( exp((1-x)(1-y)(1-z)) Uz)z = F
! 1. -(exp(-xyz) Ux)x - (exp(xyz) Uy)y  -( exp((1-x)(1-y)(1-z)) Uz)z 
!    +20sin(pi y) Ux + (sin(pi y) U)x  
!    +  sin(pi z) Uy + (sin(pi z) U)y 
!    +  sin(pi x) Uz + (sin(pi x) U)z 
!    +  1/(1+x+y+z) U                                                = F
! 2. as 1. but with Ux terms scaled by 50
! 3. as 1. but with reactive term (U) scaled by 100
! 4. as 1. but with Ux and Uy scaled by 50 and 1000, resp.

!boundary conditions: Dirichlet                                                         
!=======================================================================================
! (C) Quantum physics benchmarks                                                        
! - 3D model of Anderson localization                                                   
! 1. dx^2 Du - (4+/-rnd[L])u = F, with random numbers in [-L/2 L/2] on the diagonal     
!                                 and periodic BC                                       
!                                 (Schenk/Bollhoefer/Roemer SIAM Review 2005)           
!                                                                                       
! TODO: the right-hand sides F that follow from prescribed analytic solutions U (in   
! e.g. the A and B benchmarks) are not implemented, if we want to assess accuracy we    
! should do that, but for now we just use b=A*x for some given x.
module matpde3d_module
  implicit none

  public :: MATPDE3D_rowFunc,MATPDE3D_rhsFunc, MATPDE3D_solFunc
  public :: MATPDE3D_initDimensions, MATPDE3D_selectProblem

  real(kind=8), parameter :: pi = 4*atan(1.0_8)

  !! flags to select the test problem, encoded using 2-digit hex numbers
  integer, parameter :: PROB_A0=INT(Z'A0')
  integer, parameter :: PROB_A1=INT(Z'A1')
  integer, parameter :: PROB_A2=INT(Z'A2')
  integer, parameter :: PROB_A3=INT(Z'A3')
  integer, parameter :: PROB_A4=INT(Z'A4')
  integer, parameter :: PROB_A5=INT(Z'A5')
  integer, parameter :: PROB_A6=INT(Z'A6')
  integer, parameter :: PROB_A7=INT(Z'A7')
  integer, parameter :: PROB_A8=INT(Z'A8')
  integer, parameter :: PROB_A9=INT(Z'A9')

  integer, parameter :: PROB_B0=INT(Z'B0')
  integer, parameter :: PROB_B1=INT(Z'B1')
  integer, parameter :: PROB_B2=INT(Z'B2')
  integer, parameter :: PROB_B3=INT(Z'B3')
  integer, parameter :: PROB_B4=INT(Z'B4')
  integer, parameter :: PROB_B5=INT(Z'B5')
  integer, parameter :: PROB_B6=INT(Z'B6')
  integer, parameter :: PROB_B7=INT(Z'B7')
  integer, parameter :: PROB_B8=INT(Z'B8')
  integer, parameter :: PROB_B9=INT(Z'B9')

  integer, parameter :: PROB_C0=INT(Z'C0')
  integer, parameter :: PROB_C1=INT(Z'C1')
  integer, parameter :: PROB_C2=INT(Z'C2')
  integer, parameter :: PROB_C3=INT(Z'C3')
  integer, parameter :: PROB_C4=INT(Z'C4')
  integer, parameter :: PROB_C5=INT(Z'C5')
  integer, parameter :: PROB_C6=INT(Z'C6')
  integer, parameter :: PROB_C7=INT(Z'C7')
  integer, parameter :: PROB_C8=INT(Z'C8')
  integer, parameter :: PROB_C9=INT(Z'C9')

  ! location of boundaries in the BNDRY array
  integer, parameter :: BOTTOM=1
  integer, parameter :: TOP=2
  integer, parameter :: SOUTH=3
  integer, parameter :: NORTH=4
  integer, parameter :: WEST=5
  integer, parameter :: EAST=6


  
  !! internal switch to select one of the following test problems
  !! 1: Anderson model, -1 on off-diagonals and random numbers on diagonal
  !! 2: convection-diffusion problem with heterogenous coefficients,
  !!    
  integer :: problem
  
  integer :: BNDRY( 6 )

    ! scaling factors
    ! reaction term t(x)*u scaled by alpha,
    ! convective terms (x,y,z-direction) scaed by (beta,gamma,delta)
    real(kind=8) :: alpha,beta, gamma, delta

  integer(kind=8) :: nx, ny, nz
  
  ! grid-size dependent parameters
  real(kind=8) :: hx,hy,hz,hy2,hz2,ra,rb,raz,rbz
  
  ! for the octree ordering
  integer :: level

contains

  !! copied from the builtin/env_module.f90 because we need some
  !! random numbers for the Anderson model and the kernel lib   
  !! may not initialize the Fortran random number generator.
  subroutine init_random_seed()
    use iso_fortran_env, only: int64
#ifdef __INTEL_COMPILER
    use ifport, only: getpid
#endif
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed

  subroutine MATPDE3D_initDimensions(new_nx, new_ny, new_nz, nrows, maxnne_per_row) bind(C,name='MATPDE3D_initDimensions')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: new_nx, new_ny, new_nz
    integer(kind=G_GIDX_T), intent(out) :: nrows
    integer(kind=G_LIDX_T), intent(out) :: maxnne_per_row

    if( new_nx .ne. new_ny .or. new_nx .ne. new_nz) then
      write(*,*) 'MATPDE3D: error, nx != ny or nx != nz!'
      call exit(1)
    end if
    level = nint(log(1._8*new_nx)/log(2._8))
    if( 2**level .ne. new_nx ) then
      write(*,*) 'MATPDE3D: error, nx needs to be a power of 2!'
      call exit(1)
    end if

    nx = new_nx
    ny = new_ny
    nz = new_nz

    nrows = nx*ny*nz

    maxnne_per_row = 7

    !     set up other parameters
    HX = 1. / (DBLE(NX) + 1.)
    HY = 1. / (DBLE(NY) + 1.)
    HZ = 1. / (DBLE(NZ) + 1.)
    HY2 = HY*HY
    HZ2 = HZ*HZ
    RA = (HY/HX)**2
    RB  = HY2/HX
    RAZ = (HY/HZ)**2
    RBZ = (HY2/HZ)


  end subroutine MATPDE3D_initDimensions

  ! select problem setup (see description above).
  ! Boundary conditions are currently preset here,
  ! we might offer an additional function selectBC
  ! at some point
  !
  ! The problem is chosen by providing
  ! a hex number like e.g. 0xA3 or 0xB1 to select problem A1 or B3, respectively
  subroutine MATPDE3D_selectProblem(prob,iflag) bind(C,name='MATPDE3D_selectProblem')
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value :: prob
  integer(kind=C_INT), intent(inout) :: iflag
  
  problem = prob
  iflag=0
  ! boundary conditions are encoded like this:
  ! 0: Dirichlet
  ! 1: Neumann
  !-1: periodic
  if (problem.ge.PROB_A0 .and. PROBLEM.le.PROB_A9) then
    BNDRY(1:6)=0
  else  if (problem.ge.PROB_B0 .and. PROBLEM.le.PROB_B9) then
    BNDRY(1:6)=0
  else if (problem==PROB_C1) then
    BNDRY(1:6)=-1
    call init_random_seed()
  else
    iflag=-99
  end if
  
    ! scaling factors of convective terms (x,y,z-direction)
    alpha=1.0_8
    beta=1.0_8
    gamma=1.0_8
    delta=1.0_8
    
    if (problem==PROB_A0 .or. problem==PROB_B0) then
      alpha=0.0_8
      beta=0.0_8
      gamma=0.0_8
      delta=0.0_8
    else if (problem==PROB_A8) then
      beta=10.0_8
      gamma=10.0_8
    else if (problem==PROB_A9) then
      beta=1000.0_8
      gamma=1000.0_8
    else if (problem==PROB_B1) then
      beta=20.0_8
    else if (problem==PROB_B2) then
      beta=1000.0_8
    else if (problem==PROB_B3) then
      beta=20.0_8
      alpha=100.0_8
    else if (problem==PROB_B4) then
      beta =1000.0_8
      gamma=1000.0_8
    else if (problem .ge. PROB_C0 .and. PROBLEM .le. PROB_C9) then
      ! scaling of random numbers on diagonal
      alpha=16.5
    end if


  
  end subroutine MATPDE3D_selectProblem

  ! octree ordering
  pure function idOfCoord(coord) result(id)
    integer(kind=8), intent(in) :: coord(3)
    integer(kind=8) :: id
    integer(kind=8) :: fak2, fak8
    integer(kind=8) :: upperBound
    integer(kind=8) :: myCoord(3)
    integer :: i

    upperBound = 2**level
    myCoord = modulo(coord, upperBound)
    WHERE(myCoord==-1) myCoord=upperBound-1

    fak2 = 1
    fak8 = 1
    id = 0
    do i = 1, level
      id = id + fak8 * ( 1 * mod(myCoord(1)/fak2,2_8) &
        &              + 2 * mod(myCoord(2)/fak2,2_8) &
        &              + 4 * mod(myCoord(3)/fak2,2_8) )
      fak2 = fak2 * 2
      fak8 = fak8 * 8
    end do
  end function idOfCoord


  pure function coordOfId(id) result(coord)
    integer(kind=8), intent(in) :: id
    integer(kind=8) :: coord(3)
    integer(kind=8) :: tElem, fak(3)
    integer :: bitlevel, i

    fak(1) = 1
    fak(2) = 2
    fak(3) = 4

    bitlevel = 1
    tElem = mod(id, 8**level)

    coord = 0
    do while( (tElem / fak(1)) > 0 )
      do i = 1, 3
        coord(i) = coord(i) + bitlevel * int(mod(tElem/fak(i),2_8))
      end do
      bitlevel = bitlevel * 2
      fak = fak * 8
    end do
  end function coordOfId


  function MATPDE3D_rowFunc(row, nnz, cols, vals) result(the_result) bind(C,name='MATPDE3D_rowFunc')
    use, intrinsic :: iso_c_binding
    integer(G_GIDX_T), value :: row
    integer(G_LIDX_T), intent(inout) :: nnz
    integer(G_GIDX_T), intent(inout) :: cols(*)
    real(C_DOUBLE),    intent(inout) :: vals(*)
    integer(C_INT) :: the_result

! copied from MatrixMarket, original description:
    !
    !  Purpose
    !  =======
    !
    !  Forms the 7-point central differences operator for the elliptic
    !  PDE
    !
    !   -(P Ux)x -(Qa Uy)y  -(Qb Uz)z +  R Ux + (R U)x  +  S Uy + (S U)y + Sb Uz + (Sb U)z + T U = F
    !
    !  where P, Q(b), R, S(b) and T are the functions of x and y.  The domain
    !  is the unit cube (0,1)x(0,1)x(0,1), and the boundary condition are
    !  Dirichlet.
    !
    !  The matrix is a block tridiagonal matrix, there are NY blocks of
    !  size NX by NX on the diagonal (each block is a tridiagonal matrix),
    !  and then NY-1 blocks of size NX by NX on the sub- and super-block
    !  diagonal positions.
    !
    !  Important note: the matrix A is stored in compressed ROW format.
    !
    !  Arguments
    !  =========
    !
    !  N       (output) INTEGER
    !          The order of matrix A.  N = NX*NY
    !
    !  NX      (input) INTEGER
    !          The number of mesh points in the X-axis
    !
    !  NY      (input) INTEGER
    !          The number of mesh points in the Y-axis
    !
    !  A       (output) DOUBLE PRECISION array, dimension (NZ)
    !          The numerical values of the nonzero elements of matrix.
    !          where NZ is the number of nonzero elements in the generated
    !          matrix. Specifically, NZ = NX*NY + (NX-1)*NY*2 + NX*(NY-1)*2.
    !
    !  ROWPTR  (output) INTEGER array, dimension (N+1)
    !          The row start pointers.
    !
    !  COLIND  (output) INTEGER array, dimension (NZ)
    !          the column indices.
    !
    !  RHS     (output) DOUBLE PRECISION array, dimension (N)
    !          The right hand side b of the linear system Ax = b.
    !
    !
    !  Internal parameters
    !  ===================
    !
    !  BETA, GAMMA are the coefficients of the elliptic operator,
    !  BNDRY( 6 ) are the boundary conditions. These parameters
    !  may be changed, which will result in different matrices and spectral
    !  distribution.
    !
    !  ==================================================================
    !     ..
    !     .. Scalar variables ..
    integer(kind=8) :: ix, jy, kz, index, jd, coord(3), coord_(3)
    real(kind=8) coef
    real(kind=8) :: p12, pm12, q12, qm12, xi, rij, r1, rm1, sij, s1, sm1, yj
    real(kind=8) :: qb12, qbm12, sb1, sbij, sbm1, sbmij
    real(kind=8) :: r,zk
    
    nnz = 0

    coord = coordOfId(row)
    ix = coord(1)+1
    jy = coord(2)+1
    kz = coord(3)+1
      
    ZK = HZ*DBLE(KZ)
    YJ = HY*DBLE(JY)
    XI = HX*DBLE(IX)

    P12  = PC (XI + 0.5*HX,YJ, ZK)
    PM12 = PC (XI - 0.5*HX,YJ, ZK)
    Q12  = QC (XI,YJ + 0.5*HY, ZK)
    QM12 = QC (XI,YJ - 0.5*HY, ZK)
    QB12  = QBC (XI,YJ, ZK + 0.5*HZ)
    QBM12 = QBC (XI,YJ, ZK - 0.5*HZ)
    RIJ  = RC (BETA,XI,YJ, ZK)
    R1   = RC (BETA,XI + HX,YJ, ZK)
    RM1  = RC (BETA,XI - HX,YJ, ZK)
    SIJ  = SC (GAMMA,XI,YJ, ZK)
    S1   = SC (GAMMA,XI,YJ + HY, ZK)
    SM1  = SC (GAMMA,XI,YJ - HY, ZK)
    SBIJ  = SBC (DELTA,XI,YJ, ZK)
    SB1   = SBC (DELTA,XI,YJ, ZK + HZ)
    SBM1  = SBC (DELTA,XI,YJ, ZK - HZ)


    !           DIAGONAL.
    nnz = nnz + 1
    vals(nnz) = ra*(p12 + pm12)  +  q12 + qm12  + raz*(qb12+qbm12) + hy2*tc(alpha,xi,yj,zk)
    cols(nnz) = idOfCoord(coord) !index
    jd = nnz
    
    if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! add random diagonal term
      call random_number(r)
      vals(nnz)=vals(nnz)+alpha*(r-0.5_8)
    end if

    ! LOWEST TWO BANDS
    if (kz.ne.1 .or. BNDRY(BOTTOM)==-1) then
      nnz = nnz + 1
      vals(nnz) = - ( raz*qbm12 + 0.5*rbz*(sbij+sbm1) )
      coord_ = coord-(/0,0,1/)
      cols(nnz) = idOfCoord(coord_)
    end if

    if (jy.ne.1 .or. BNDRY(SOUTH)==-1) then
      nnz = nnz + 1
      vals(nnz) = - ( qm12 + 0.5*hy*(sij+sm1) )
      coord_ = coord-(/0,1,0/)
      cols(nnz) = idOfCoord(coord_) 
    end if

    !           SUB-DIAGONAL.
    if (ix.ne.1 .or. BNDRY(WEST)==-1) then
      nnz = nnz + 1
      vals(nnz) = - ( ra*pm12 + 0.5*rb*(rij+rm1) )
      coord_ = coord-(/1,0,0/)
      cols(nnz) = idOfCoord(coord_) !index - 1
    end if

    !           SUPER-DIAGONAL.
    if (ix.ne.nx .or. BNDRY(EAST)==-1)  then
      nnz = nnz + 1
      vals(nnz) = -ra*p12 + 0.5*rb*(rij+r1)
      coord_ = coord+(/1,0,0/)
      cols(nnz) = idOfCoord(coord_) !index + 1
    end if

    !           HIGHEST BANDS.
    if (jy.ne.ny .or. BNDRY(NORTH)==-1) then
      nnz = nnz + 1
      vals(nnz) = -q12 + 0.5*hy*(sij+s1)
      coord_ = coord+(/0,1,0/)
      cols(nnz) = idOfCoord(coord_) !index + nx
    end if

    if (kz.ne.nz .or. BNDRY(top)==-1) then
      nnz = nnz + 1
      vals(nnz) = -raz*qb12 + 0.5*rbz*(sbij+sb1)
      coord_ = coord+(/0,0,1/)
      cols(nnz) = idOfCoord(coord_) !index + nx
    end if

    !           BOUNDARY CONDITIONS.
    if (kz.eq.1) then
      COEF = - ( RAZ*QBM12 + 0.5*RBZ*(SBIJ+SBM1) )
      IF (BNDRY(BOTTOM).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (jy.eq.1) then
      COEF = - ( QM12 + 0.5*HY*(SIJ+SM1) )
      IF (BNDRY(SOUTH).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (ix.eq.1) then
      COEF = - ( RA*PM12 + 0.5*RB*(RIJ+RM1) )
      IF (BNDRY(WEST).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (IX.eq.NX) then
      COEF = -RA*P12 + 0.5*RB*(RIJ+R1)
      IF (BNDRY(EAST).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (JY.eq.NY) then
      COEF = -Q12 + 0.5*HY*(SIJ+S1)
      if (BNDRY(NORTH).eq.1) vals(jd) = vals(jd) + coef
    end if
    IF (KZ.eq.NZ) then
      COEF = -RAZ*QB12 + 0.5*RBZ*(SBIJ+SB1)
      if (BNDRY(TOP).eq.1) vals(jd) = vals(jd) + coef
    end if

    the_result = 0
    
  end function MATPDE3D_rowFunc


  !! construct the RHS for a linear system AU=F with A
  !! as provided by MATPDE3D_rowFunc and U as provided
  !! by MATPDE3D_solFunc. The same initialization procedure
  !! as for MATPDE3D_rowFunc is required, i.e. 
  !! first initDimensions and then selectProblem.
  function MATPDE3D_rhsFunc(row, rhsVal,col) result(the_result) bind(C,name='MATPDE3D_rhsFunc')
    use, intrinsic :: iso_c_binding
    integer(G_GIDX_T), value :: row
    integer(G_LIDX_T), value :: col
    real(C_DOUBLE),    intent(inout) :: rhsVal
    integer(C_INT) :: the_result

    !     .. Scalar variables ..
    integer(kind=8) :: ix, jy, kz, index, jd, coord(3)
    real(kind=8) coef
    real(kind=8) :: xi,yj,zk

    real(kind=8) :: p12, pm12, q12, qm12, rij, r1, rm1, sij, s1, sm1
    real(kind=8) :: qb12, qbm12, sb1, sbij, sbm1, sbmij
    
    coord = coordOfId(row)
    ix = coord(1)+1
    jy = coord(2)+1
    kz = coord(3)+1
      
    ZK = HZ*DBLE(KZ)
    YJ = HY*DBLE(JY)
    XI = HX*DBLE(IX)

    if (ix.eq.1.or.ix.eq.nx .or. &
        jy.eq.1.or.jy.eq.ny .or. &
        kz.eq.1.or.kz.eq.nz ) then
        
    P12  = PC (XI + 0.5*HX,YJ, ZK)
    PM12 = PC (XI - 0.5*HX,YJ, ZK)
    Q12  = QC (XI,YJ + 0.5*HY, ZK)
    QM12 = QC (XI,YJ - 0.5*HY, ZK)
    QB12  = QBC (XI,YJ, ZK + 0.5*HZ)
    QBM12 = QBC (XI,YJ, ZK - 0.5*HZ)
    RIJ  = RC (BETA,XI,YJ, ZK)
    R1   = RC (BETA,XI + HX,YJ, ZK)
    RM1  = RC (BETA,XI - HX,YJ, ZK)
    SIJ  = SC (GAMMA,XI,YJ, ZK)
    S1   = SC (GAMMA,XI,YJ + HY, ZK)
    SM1  = SC (GAMMA,XI,YJ - HY, ZK)
    SBIJ  = SBC (DELTA,XI,YJ, ZK)
    SB1   = SBC (DELTA,XI,YJ, ZK + HZ)
    SBM1  = SBC (DELTA,XI,YJ, ZK - HZ)
    
    end if

    rhsVal = hy2*f(xi,yj,zk,col)
    
    !           BOUNDARY CONDITIONS.
    if (kz.eq.1 .AND. BNDRY(BOTTOM).EQ.0) THEN
      COEF = - ( RAZ*QBM12 + 0.5*RBZ*(SBIJ+SBM1) )
      rhsVal = rhsVal - COEF*U(XI,YJ-HY,ZK,col)
    end if
    if (jy.eq.1 .AND. BNDRY(SOUTH).EQ.0) THEN
      COEF = - ( QM12 + 0.5*HY*(SIJ+SM1) )
      rhsVal = rhsVal - COEF*U(XI,YJ-HY,ZK,col)
    end if
    if (ix.eq.1 .AND. BNDRY(WEST).EQ.0) THEN
      COEF = - ( RA*PM12 + 0.5*RB*(RIJ+RM1) )
      rhsVal = rhsVal - COEF*U(XI-HX,YJ,ZK,col)
    end if
    IF (IX.eq.NX .AND. BNDRY(EAST).EQ.0) then
      COEF = -RA*P12 + 0.5*RB*(RIJ+R1)
      rhsVal = rhsVal - COEF*U(XI+HX,YJ,ZK,col)
    end if
    IF (JY.eq.NY .AND. BNDRY(NORTH).eq.0) then
      COEF = -Q12 + 0.5*HY*(SIJ+S1)
      rhsVal = rhsVal - COEF*U(XI,YJ+HY,ZK,col)
    end if
    IF (KZ.eq.NZ .AND. BNDRY(TOP).eq.0) then
      COEF = -RAZ*QB12 + 0.5*RBZ*(SBIJ+SB1)
      rhsVal = rhsVal - COEF*U(XI,YJ+HY,ZK,col)
    end if

    the_result = 0
    
  end function MATPDE3D_rhsFunc

  function MATPDE3D_solFunc(row, solVal,col) result(the_result) bind(C,name='MATPDE3D_solFunc')
    use, intrinsic :: iso_c_binding
    integer(G_GIDX_T), value :: row
    integer(G_LIDX_T), value :: col
    real(C_DOUBLE),    intent(inout) :: solVal
    integer(C_INT) :: the_result

    !     .. Scalar variables ..
    integer(kind=8) :: ix, jy, kz, index, jd, coord(3)
    real(kind=8) :: xi,yj,zk
    
    coord = coordOfId(row)
    ix = coord(1)+1
    jy = coord(2)+1
    kz = coord(3)+1
      
    ZK = HZ*DBLE(KZ)
    YJ = HY*DBLE(JY)
    XI = HX*DBLE(IX)
    
    solVal=U(XI,YJ,ZK,col)
    the_result=0
    
  end function MATPDE3D_solFunc

  pure function pc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: pc
    
    if (problem .ge. PROB_A0 .and. problem .le. PROB_A9) then
      ! Gordon problems A 1-9
      pc = 1.0_8
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      pc = exp(-x*y*z)
    else if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      pc = 1.0_8
    else
      ! default
      pc = 1.0_8
    end if

  end function pc

  ! for constructing the RHS for a given U
  pure function dpdx(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: dpdx
    
    if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      dpdx = -y*z*exp(-x*y*z)
    else
      ! default: constant p
      dpdx = 0.0_8
    end if

  end function dpdx

  pure function qc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: qc

    if (problem .ge. PROB_A0 .and. problem .le. PROB_A9) then
      ! Gordon problems A1-9
      qc = 1.0_8
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      qc = exp(x*y*z)
    else if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      qc = 1.0_8
    else
      ! default
      qc = 1.0_8
    end if

  end function qc

  pure function dqdy(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: dqdy

    if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      dqdy = x*z*exp(x*y*z)
    else
      dqdy = 0.0_8
    end if

  end function dqdy

  pure function qbc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: qbc

    if (problem .ge. PROB_A0 .and. problem .le. PROB_A9) then
      ! Gordon problems A 1-9
      qbc = 1.0_8
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      qbc = exp((1.0_8-x)*(1.0_8-y)*(1.0_8-z))
    else if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      qbc = 1.0_8
    else
      ! default
      qbc = 1.0_8
    end if

  end function qbc

  pure function dqbdz(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: dqbdz

    if (problem .ge. PROB_A0 .and. problem .le. PROB_A9) then
      ! Gordon problems A 1-9
      dqbdz = 1.0_8
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      dqbdz = -((1.0_8-x)*(1.0_8-y))*exp((1.0_8-x)*(1.0_8-y)*(1.0_8-z))
    else
      ! default
      dqbdz = 0.0_8
    end if

  end function dqbdz

  pure function rc(beta, x,y, z)
    real(kind=8), intent(in) :: beta, x, y, z
    real(kind=8) :: rc
    
    if (problem == PROB_A1) then
      rc = 500.0_8
    else if (problem == PROB_A2) then
      rc = 500.0_8*exp(x*y*z)
    else if (problem == PROB_A3) then
      rc = 50.0_8*x
    else if (problem == PROB_A4) then
      rc = -0.5e5*x*x
    else if (problem == PROB_A5) then
      rc = -500.0_8*(1.0_8+x*x)
    else if (problem == PROB_A6) then
      rc = -500.0_8*(1.0_8-2.0_8*x)
    else if (problem == PROB_A7) then
      rc = -500.0_8*x*x
    else if (problem == PROB_A8 .or. problem == PROB_A9) then
      rc = -0.5_8*exp(x*y)
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      rc = sin(pi*y)
    else if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      rc = 0.0_8
    else
      ! default
      rc = 0.0_8
    end if    
    rc=rc*beta
  end function rc

  pure function drdx(beta, x,y, z)
    real(kind=8), intent(in) :: beta, x, y, z
    real(kind=8) :: drdx
    
    if (problem == PROB_A2) then
      drdx = 500.0_8*y*z*exp(x*y*z)
    else if (problem == PROB_A3) then
      drdx = 50.0_8
    else if (problem == PROB_A4) then
      drdx = -1.0e5*x
    else if (problem == PROB_A5) then
      drdx = -1000.0_8*x
    else if (problem == PROB_A6) then
      drdx = 1000.0_8
    else if (problem == PROB_A7) then
      drdx = -1000.0_8*x
    else if (problem == PROB_A8 .or. problem == PROB_A9) then
      drdx = -0.5_8*y*exp(x*y)
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      drdx = sin(pi*y)
    else if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      drdx = 0.0_8
    else
      ! default
      drdx = 0.0_8
    end if    
    drdx=drdx*beta
  end function drdx

  pure function sc(gamma, x,y, z)
    real(kind=8), intent(in) :: gamma, x, y, z
    real(kind=8) :: sc

    if (problem == PROB_A1) then
      sc = 0.0_8
    else if (problem == PROB_A2) then
      sc = 500.0_8*exp(x*y*z)
    else if (problem == PROB_A3) then
      sc = -0.5_8*y
    else if (problem == PROB_A4) then
      sc = -0.5e5*x*x
    else if (problem == PROB_A5) then
      sc = 50.0_8
    else if (problem == PROB_A6) then
      sc = -500.0_8*(1.0_8-2.0_8*y)
    else if (problem == PROB_A7) then
      sc = 0.0_8
    else if (problem==PROB_A8 .or. problem==PROB_A9) then
      sc = -0.5_8*exp(-x*y)
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      sc = sin(pi*z)
    else
      ! default
      sc = 0.0_8
    end if
    
    sc=sc*gamma

  end function sc

  pure function dsdy(gamma, x,y, z)
    real(kind=8), intent(in) :: gamma, x, y, z
    real(kind=8) :: dsdy

    if (problem == PROB_A2) then
      dsdy = 500.0_8*x*z*exp(x*y*z)
    else if (problem == PROB_A3) then
      dsdy = -0.5_8
    else if (problem == PROB_A6) then
      dsdy = 1000.0_8
    else if (problem==PROB_A8 .or. problem==PROB_A9) then
      dsdy = -0.5_8*x*z*exp(-x*y)
    else
      ! default
      dsdy = 0.0_8
    end if
    
    dsdy=dsdy*gamma

  end function dsdy

  pure function sbc(delta, x,y,z)
    real(kind=8), intent(in) :: delta, x, y, z
    real(kind=8) :: sbc

    if (problem == PROB_A1) then
      sbc = 0.0_8
    else if (problem == PROB_A2) then
      sbc = -500.0_8*exp(x*y*z)
    else if (problem == PROB_A3) then
      sbc = +0.5_8*z
    else if (problem == PROB_A4) then
      sbc = -0.5e5*x*x
    else if (problem == PROB_A5) then
      sbc = 50.0_8
    else if (problem == PROB_A6) then
      sbc = -500.0_8*(1.0_8-2.0_8*z)
    else if (problem == PROB_A7) then
      sbc = 0.0_8
    else if (problem == PROB_A8) then
      sbc = 0.0_8
    else if (problem == PROB_A9) then
      sbc = 0.0_8
    else if (problem .ge. PROB_B0 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      sbc = sin(pi*x)
    else if (problem .ge. PROB_C0 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      sbc = 0.0_8
    else
      ! default
      sbc = 0.0_8
    end if
    
    sbc=sbc*delta

  end function sbc

  pure function dsbdz(delta, x,y,z)
    real(kind=8), intent(in) :: delta, x, y, z
    real(kind=8) :: dsbdz

    if (problem == PROB_A2) then
      dsbdz = -500.0_8*x*y*exp(x*y*z)
    else if (problem == PROB_A3) then
      dsbdz = +0.5_8
    else if (problem == PROB_A6) then
      dsbdz = -1000.0_8
    else
      ! default
      dsbdz = 0.0_8
    end if
    
    dsbdz=dsbdz*delta

  end function dsbdz

  pure function tc(alpha,x,y,z)
    real(kind=8), intent(in) :: alpha
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: tc
  !term in front of u.
  ! Note that some terms appear here because
  ! in the Gordon paper they use e.g. ru_x
  ! whereas in matpde we assume (r u)_x + r u_x
  if (problem==PROB_A2) then
    tc = (y*z+x*z-x*y)*500
  else if (problem==PROB_A3) then
    tc = 50.0_8+100.0_8*(x+y+z)/(x*y*z)
  else if (problem==PROB_A4) then
    tc = 1.0e5
  else if (problem==PROB_A5) then
    tc=-1000.0_8*x
  else if (problem==PROB_A6) then
    tc=3000.0_8
  else if (problem==PROB_A7) then
    tc = 0.0_8
  else if (problem==PROB_A8) then
    tc = -5.0_8*(exp(x*y)+exp(-x*y))
  else if (problem==PROB_A9) then
    tc = -500.0_8*(exp(x*y)+exp(-x*y))
  else if (problem == PROB_B1) then
    tc = 1. / (1.+x+y+z)
  else if (problem == PROB_C1) then
    tc = -6.0_8/(hy2*alpha)
  else
    tc = 0.0_8
  end if

  tc=tc*alpha

  end function tc

  ! Analytical solution we prescribe
  ! the derivatives were calculated using the matlab lines
  !
  !syms x y z k pi
  !u    = x*exp(x*y*z)*sin(k*pi*x)*sin(k*pi*y)*sin(k*pi*z);
  !ux = simplify(diff(u,x)); uy = simplify(diff(u,y)); uz=simplify(diff(u,z));
  !uxx= simplify(diff(ux,x)); uyy=simplify(diff(uy,y)); uzz= simplify(diff(uz,z));
  pure function u(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: u
  
    if (problem==PROB_A1) then
      u = x * x*y*z*(1.-x)*(1.-y)*(1.-z)
    else if (problem==PROB_A2) then
      u = x+y+z
    else
      u = x * exp(x*y*z) * sin(k*pi*x) * sin(k*pi*y) * sin(k*pi*z)
    end if
  
  end function u

  pure function ux(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: ux, cx,cy,cz,sx,sy,sz,exyz
    
    if (problem==PROB_A1) then
      if (x==0. .or. x==1.) then
        ux=0.
      else
        ux = u(x,y,z,k)/x - u(x,y,z,k)/(1.-x)
      end if
    else if (problem==PROB_A2) then
      ux=1.
    else
      sx = sin(k*pi*x)
      cx = cos(k*pi*x)
      sy = sin(k*pi*y)
      cy = cos(k*pi*y)
      sz = sin(k*pi*z)
      cz = cos(k*pi*z)
      exyz=exp(x*y*z)
      ux=sy*sz*exyz*(sx + pi*k*x*cx + x*y*z*sx)
    end if
  
  end function ux

  pure function uy(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: uy, cx,cy,cz,sx,sy,sz,exyz

    if (problem==PROB_A1) then
      if (y==0. .or. y==1.) then
        uy=0.
      else
        uy = u(x,y,z,k)/y - u(x,y,z,k)/(1.-y)
      end if
    else if (problem==PROB_A2) then
      uy=1.
    else
      sx = sin(k*pi*x)
      cx = cos(k*pi*x)
      sy = sin(k*pi*y)
      cy = cos(k*pi*y)
      sz = sin(k*pi*z)
      cz = cos(k*pi*z)
      exyz=exp(x*y*z)

      uy=x*sx*sz*exyz*(pi*k*cy + x*z*sy)
    end if
  end function uy

  pure function uz(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: uz, cx,cy,cz,sx,sy,sz,exyz


    if (problem==PROB_A1) then
      if (z==0. .or. z==1.) then
        uz=0.
      else
        uz = u(x,y,z,k)/z - u(x,y,z,k)/(1.-z)
      end if
    else if (problem==PROB_A2) then
      uz=1.
    else
      sx = sin(k*pi*x)
      cx = cos(k*pi*x)
      sy = sin(k*pi*y)
      cy = cos(k*pi*y)
      sz = sin(k*pi*z)
      cz = cos(k*pi*z)
      exyz=exp(x*y*z)
      uz=x*sx*sy*exyz*(pi*k*cz + x*y*sz)
    end if
  
  end function uz

  pure function uxx(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: uxx, cx,cy,cz,sx,sy,sz,exyz


  if (problem==PROB_A1.or.PROBLEM==PROB_A2) then
    uxx=0.0
  else
    sx = sin(k*pi*x)
    cx = cos(k*pi*x)
    sy = sin(k*pi*y)
    cy = cos(k*pi*y)
    sz = sin(k*pi*z)
    cz = cos(k*pi*z)
    exyz=exp(x*y*z)
    
uxx=sy*sz*exyz*(2.0_8*pi*k*cx + 2*y*z*sx - (pi*k)**2.0*x*sx + x*y*y*z*z*sx + 2.0*pi*k*x*y*z*cx)
end if
  
  end function uxx

  pure function uyy(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: uyy, cx,cy,cz,sx,sy,sz,exyz

  if (problem==PROB_A1.or.PROBLEM==PROB_A2) then
    uyy=0.0
  else
    sx = sin(k*pi*x)
    cx = cos(k*pi*x)
    sy = sin(k*pi*y)
    cy = cos(k*pi*y)
    sz = sin(k*pi*z)
    cz = cos(k*pi*z)
    exyz=exp(x*y*z)

uyy=x*sx*sz*exyz*(x*x*z*z*sy - (pi*k)**2*sy + 2.0*pi*k*x*z*cy)
end if 
 
  end function uyy

  pure function uzz(x,y,z,k)
  real(kind=8), intent(in) :: x, y, z
  integer(kind=4), intent(in) :: k
  real(kind=8) :: uzz, cx,cy,cz,sx,sy,sz,exyz

  if (problem==PROB_A1.or.PROBLEM==PROB_A2) then
    uzz=0.0
  else
    sx = sin(k*pi*x)
    cx = cos(k*pi*x)
    sy = sin(k*pi*y)
    cy = cos(k*pi*y)
    sz = sin(k*pi*z)
    cz = cos(k*pi*z)
    exyz=exp(x*y*z)
  
uzz=x*sx*sy*exyz*(x*x*y*y*sz - (pi*k)**2*sz + 2.0*pi*k*x*y*cz)
end if

  end function uzz

  pure function f(x,y,z,k)
    real(kind=8), intent(in) :: x, y, z
    integer, intent(in) :: k
    real(kind=8) :: f

    real(kind=8) :: a, ax, axx, ay, ayy, az, azz
    real(kind=8) p, px, q, qy, qb,qbz, r, rx,s,sy,sb,sbz,t
    
    f = 0.0_8

    ! evaluate analytic solution and its derivatives
    a=u(x,y,z,k)
    ax=ux(x,y,z,k)
    ay=uy(x,y,z,k)
    az=uz(x,y,z,k)
    axx=uxx(x,y,z,k)
    ayy=uyy(x,y,z,k)
    azz=uzz(x,y,z,k)
    
    ! remaining terms to evaluate F
    p   =  pc(x,y,z)
    px  =  dpdx(x,y,z)
    q   =  qc(x,y,z)
    qy  =  dqdy(x,y,z)
    qb  =  qbc(x,y,z)
    qbz =  dqbdz(x,y,z)
    r   =  rc(beta,x,y,z)
    rx  = drdx(beta,x,y,z)
    s   =  sc(gamma,x,y,z)
    sy  = dsdy(gamma,x,y,z)
    sb  = sbc(delta,x,y,z)
    sbz = dsbdz(delta,x,y,z)
    t   =  tc(alpha,x,y,z)

    f = - (p*axx + px*ax + q*ayy + qy*ay +qb*az + qbz*azz) + &
        +2.0*(r*ax+s*ay+sb*az) + (rx+sy+sbz+t)*a
    

  end function f

end module matpde3d_module
