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
! 1. -(exp(-xyz) Ux)x - (exp(xyz) Uy)y  -( exp((1-x)(1-y)(1-z)) Uz)z 
!    +  sin(pi y) Ux + (sin(pi y) U)x  
!    +  sin(pi z) Uy + (sin(pi z) U)y 
!    +  sin(pi x) Uz + (sin(pi x) U)z 
!    +  1/(1+x+y+z) U = F

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

  public :: MATPDE3D_rowFunc
  public :: MATPDE3D_initDimensions, MATPDE3D_selectProblem

  real(kind=8), parameter :: pi = 4*atan(1.0_8)

  !! flags to select the test problem, encoded using 2-digit hex numbers
  integer, parameter :: PROB_A1=INT(Z'A1')
  integer, parameter :: PROB_A2=INT(Z'A2')
  integer, parameter :: PROB_A3=INT(Z'A3')
  integer, parameter :: PROB_A4=INT(Z'A4')
  integer, parameter :: PROB_A5=INT(Z'A5')
  integer, parameter :: PROB_A6=INT(Z'A6')
  integer, parameter :: PROB_A7=INT(Z'A7')
  integer, parameter :: PROB_A8=INT(Z'A8')
  integer, parameter :: PROB_A9=INT(Z'A9')

  integer, parameter :: PROB_B1=INT(Z'B1')
  integer, parameter :: PROB_B2=INT(Z'B2')
  integer, parameter :: PROB_B3=INT(Z'B3')
  integer, parameter :: PROB_B4=INT(Z'B4')
  integer, parameter :: PROB_B5=INT(Z'B5')
  integer, parameter :: PROB_B6=INT(Z'B6')
  integer, parameter :: PROB_B7=INT(Z'B7')
  integer, parameter :: PROB_B8=INT(Z'B8')
  integer, parameter :: PROB_B9=INT(Z'B9')

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


  integer(kind=8) :: nx, ny, nz
  
  ! grid-size dependent parameters
  real(kind=8) :: hx,hy,hz,hy2,hz2,ra,rb,raz,rbz
  
  ! for the octree ordering
  integer :: level

contains

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

  ! select problem setup (see description above) and
  ! boundary conditions. The problem is chosen by providing
  ! a hex number like e.g. 0xA3 or 0xB1 to select problem A1 or B3, respectively,
  ! and the boundary conditions are given as an integer array for the various boundaries
  ! in the order [Z0 Z1 Y0 Y1 X0 X1] (e.g. X0 being the boundary at x=0).
  ! The values in the array indicate the type of boundary condition:
  !  0: homogenous Dirichlet
  !  1: Dirichlet with a value taken the analytic solution
  ! -1: periodic boundary conditions
  subroutine MATPDE3D_selectProblem(prob, bc) bind(C,name='MATPDE3D_selectProblem')
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), intent(in) :: prob
  integer(kind=C_INT), intent(in) :: bc(6)
  
  problem = prob
  BNDRY(1:6)=bc(1:6)
  
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
    !
    !     .. Parameters ..
    real(kind=8), parameter :: beta = 20., gamma = 0., gammab=0.
    !     ..
    !     .. Scalar variables ..
    integer(kind=8) :: ix, jy, kz, index, jd, coord(3), coord_(3)
    real(kind=8) :: hx, hy, hz, hy2, hz2, ra, rb, raz, rbz, coef
    real(kind=8) :: p12, pm12, q12, qm12, xi, rij, r1, rm1, sij, s1, sm1, yj
    real(kind=8) :: qb12, qbm12, sb1, sbij, sbm1, sbmij
    real(kind=8) :: zk

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
    SBIJ  = SBC (GAMMAB,XI,YJ, ZK)
    SB1   = SBC (GAMMAB,XI,YJ, ZK + HZ)
    SBM1  = SBC (GAMMAB,XI,YJ, ZK - HZ)


    !           DIAGONAL.
    nnz = nnz + 1
    vals(nnz) = ra*(p12 + pm12)  +  q12 + qm12  + raz*(qb12+qbm12) + hy2*tc(xi,yj,zk)
    cols(nnz) = idOfCoord(coord) !index
    jd = nnz

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
      !      IF (BNDRY(1).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ-HY)
      IF (BNDRY(BOTTOM).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (jy.eq.1) then
      COEF = - ( QM12 + 0.5*HY*(SIJ+SM1) )
      !      IF (BNDRY(1).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ-HY)
      IF (BNDRY(SOUTH).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (ix.eq.1) then
      COEF = - ( RA*PM12 + 0.5*RB*(RIJ+RM1) )
      !      IF (BNDRY(2).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI-HX,YJ)
      IF (BNDRY(WEST).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (IX.eq.NX) then
      COEF = -RA*P12 + 0.5*RB*(RIJ+R1)
      !      IF (BNDRY(3).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI+HX,YJ)
      IF (BNDRY(EAST).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (JY.eq.NY) then
      COEF = -Q12 + 0.5*HY*(SIJ+S1)
      !      IF (BNDRY(4).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ+HY)
      if (BNDRY(NORTH).eq.1) vals(jd) = vals(jd) + coef
    end if
    IF (KZ.eq.NZ) then
      COEF = -RAZ*QB12 + 0.5*RBZ*(SBIJ+SB1)
      !      IF (BNDRY(4).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ+HY)
      if (BNDRY(TOP).eq.1) vals(jd) = vals(jd) + coef
    end if

    the_result = 0
    
  end function MATPDE3D_rowFunc


  pure function pc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: pc
    
    if (problem .ge. PROB_A1 .and. problem .le. PROB_A9) then
      ! Gordon problems A 1-9
      pc = 1.0_8
    else if (problem .ge. PROB_B1 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      pc = exp(-x*y*z)
    else if (problem .ge. PROB_C1 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      pc = HX*HX
    else
      ! default
      pc = 1.0_8
    end if

  end function pc

  pure function qc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: qc

    if (problem .ge. PROB_A1 .and. problem .le. PROB_A9) then
      ! Gordon problems A1-9
      qc = 1.0_8
    else if (problem .ge. PROB_B1 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      qc = exp(x*y*z)
    else if (problem .ge. PROB_C1 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      qc = HY*HY
    else
      ! default
      qc = 1.0_8
    end if

  end function qc

  pure function qbc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: qbc

    if (problem .ge. PROB_A1 .and. problem .le. PROB_A9) then
      ! Gordon problems A 1-9
      qbc = 1.0_8
    else if (problem .ge. PROB_B1 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      qbc = exp((1.0_8-x)*(1.0_8-y)*(1.0_8-z))
    else if (problem .ge. PROB_C1 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      qbc = HZ*HZ
    else
      ! default
      qbc = 1.0_8
    end if

  end function qbc

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
    else if (problem == PROB_A8) then
      rc = -5.0_8*exp(x*y)
    else if (problem == PROB_A9) then
      rc = -500.0_8*exp(x*y)
    else if (problem .ge. PROB_B1 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      rc = sin(pi*y)
    else if (problem .ge. PROB_C1 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      rc = 0.0_8
    else
      ! default
      rc = 0.0_8
    end if    
  end function rc

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
    else if (problem == PROB_A8) then
      sc = -5.0_8*exp(-x*y)
    else if (problem == PROB_A9) then
      sc = -500.0_8*exp(-x*y)
    else if (problem .ge. PROB_B1 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      sc = sin(pi*z)
    else if (problem .ge. PROB_C1 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      sc = 0.0_8
    else
      ! default
      sc = 0.0_8
    end if    

  end function sc

  pure function sbc(gamma, x,y,z)
    real(kind=8), intent(in) :: gamma, x, y, z
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
    else if (problem .ge. PROB_B1 .and. problem .le. PROB_B9) then
      ! varying coefficient problems (B)
      sbc = sin(pi*x)
    else if (problem .ge. PROB_C1 .and. problem .le. PROB_C9) then
      ! QM test cases with constant -1 in off-diagonals
      sbc = 0.0_8
    else
      ! default
      sbc = 0.0_8
    end if

  end function sbc

  ! because of the random_number for Anderson localization, this one can't be pure
  function tc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: tc
    real(kind=8) :: r
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
    call random_number(r)
    tc = -6.0_8 + 16.5_8*(r-0.5_8)
  else
    tc = 0.0_8
  end if

  end function tc

  ! analytical solution we prescribe
  pure function u(x,y,z)
  real(kind=8), intent(in) :: x, y, z
  real(kind=8) :: u
    u = x * exp(x*y*z) * sin(pi*x) * sin(pi*y) * sin(pi*z)
  
  end function u

  pure function f(beta, gamma, x,y,z)
    real(kind=8), intent(in) :: beta, gamma, x, y, z
    real(kind=8) :: f

    real(kind=8) pxy, pxxy, qxy, qyxy, cxy, rxy, rxxy, sxy, syxy, exy, sx, cx, sy, cy, a, ax, axx, ay, ayy
    
    f = 0.0_8

!    pxy  = pc(x,y,z)
!    pxxy = -y * exp(-x*y)
!    qxy  = qc(x,y,z)
!    qyxy = x * exp(x*y)
!    cxy  = 1. / ( 1. + x + y)
!    rxy  = rc(beta,x,y,z)
!    rxxy = beta
!    sxy  = sc(gamma,x,y,z)
!    syxy = gamma
!
!    exy = exp(x*y)
!    sx  = sin(pi*x)
!    cx  = cos(pi*x)
!    sy  = sin(pi*y)
!    cy  = cos(pi*y)
!
!    a   = x * exy * sx * sy0
!    ax  = exy * (pi * x * cx + (1. + x * y) * sx) * sy
!    axx = exy * (pi * (-pi * x * sx + (x * y + 2.) * cx) + y * sx) * sy + y * ax
!    ay  = exy * (pi * cy + x * sy) * x * sx
!    ayy = x * (exy * pi * (-pi * sy + x * cy) * sx + ay)
!
!    f = - (pxy*axx + qxy*ayy)  + (2.*rxy - pxxy) * ax  +  (2.*sxy - qyxy) * ay  + (rxxy + syxy + cxy) * a

  end function f

end module matpde3d_module
