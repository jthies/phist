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
module matpde3d_module
  implicit none

  public :: MATPDE3D_rowFunc
  public :: MATPDE3D_initDimensions

  real(kind=8), parameter :: pi = 4*atan(1.0_8)

  integer(kind=8) :: nx, ny, nz
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

  end subroutine MATPDE3D_initDimensions


  ! TODO octree ordering
  pure function idOfCoord(coord) result(id)
    integer(kind=8), intent(in) :: coord(3)
    integer(kind=8) :: id
    integer(kind=8) :: fak2, fak8
    integer(kind=8) :: upperBound
    integer(kind=8) :: myCoord(3)
    integer :: i

    upperBound = 2**level
    myCoord = modulo(coord, upperBound)

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


  function MATPDE3D_rowFunc(row, nnz, cols, vals) result(the_result) bind(C, name='MATPDE3D_rowFunc')
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
    !   -(P Ux)x -(Qa Uy)y  -(Qb Uz)z +  R Ux + (R U)x  +  S Uy + (S U)y + Sb Uz + (Sb U)z +  T U = F
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
    !  BNDRY( 4 ) are the boundary Dirichlet conditions. These parameters
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
    real(kind=8) :: bndry( 6 )
    real(kind=8) :: hx, hy, hz, hy2, hz2, ra, rb, raz, rbz, coef
    real(kind=8) :: p12, pm12, q12, qm12, xi, rij, r1, rm1, sij, s1, sm1, yj
    real(kind=8) :: qb12, qbm12, sb1, sbij, sbm1, sbmij
    real(kind=8) :: zk

    !     set up initial boundary conditions
    BNDRY( 1 ) = 0.
    BNDRY( 2 ) = 0.
    BNDRY( 3 ) = 0.
    BNDRY( 4 ) = 0.
    BNDRY( 5 ) = 0.
    BNDRY( 6 ) = 0.

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

    nnz = 0
    ! do jy = 1, ny
    !   do ix = 1, nx
    coord = coordOfId(row)
    ix = coord(1)+1
    jy = coord(2)+1
    kz = coord(3)+1
    !jy = row/nx+1
    !ix = mod(row,nx)+1
      
    ZK = HZ*DBLE(KZ)
    YJ = HY*DBLE(JY)
    XI = HX*DBLE(IX)
    !INDEX = IX + (JY-1)*NX
    !    RHS(INDEX) = HY2*F(BETA,GAMMA,XI,YJ)
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
    if (kz.ne.1) then
      nnz = nnz + 1
      vals(nnz) = - ( raz*qbm12 + 0.5*rbz*(sbij+sbm1) )
      coord_ = coord-(/0,0,1/)
      cols(nnz) = idOfCoord(coord_) 
    end if

    if (jy.ne.1) then
      nnz = nnz + 1
      vals(nnz) = - ( qm12 + 0.5*hy*(sij+sm1) )
      coord_ = coord-(/0,1,0/)
      cols(nnz) = idOfCoord(coord_) 
    end if

    !           SUB-DIAGONAL.
    if (ix.ne.1) then
      nnz = nnz + 1
      vals(nnz) = - ( ra*pm12 + 0.5*rb*(rij+rm1) )
      coord_ = coord-(/1,0,0/)
      cols(nnz) = idOfCoord(coord_) !index - 1
    end if

    !           SUPER-DIAGONAL.
    if (ix.ne.nx)  then
      nnz = nnz + 1
      vals(nnz) = -ra*p12 + 0.5*rb*(rij+r1)
      coord_ = coord+(/1,0,0/)
      cols(nnz) = idOfCoord(coord_) !index + 1
    end if

    !           HIGHEST BANDS.
    if (jy.ne.ny) then
      nnz = nnz + 1
      vals(nnz) = -q12 + 0.5*hy*(sij+s1)
      coord_ = coord+(/0,1,0/)
      cols(nnz) = idOfCoord(coord_) !index + nx
    end if

    if (kz.ne.nz) then
      nnz = nnz + 1
      vals(nnz) = -raz*qb12 + 0.5*rbz*(sbij+sb1)
      coord_ = coord+(/0,0,1/)
      cols(nnz) = idOfCoord(coord_) !index + nx
    end if

    !           BOUNDARY CONDITIONS.
    if (kz.eq.1) then
      COEF = - ( RAZ*QBM12 + 0.5*RBZ*(SBIJ+SBM1) )
      !      IF (BNDRY(1).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ-HY)
      IF (BNDRY(1).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (jy.eq.1) then
      COEF = - ( QM12 + 0.5*HY*(SIJ+SM1) )
      !      IF (BNDRY(1).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ-HY)
      IF (BNDRY(2).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (ix.eq.1) then
      COEF = - ( RA*PM12 + 0.5*RB*(RIJ+RM1) )
      !      IF (BNDRY(2).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI-HX,YJ)
      IF (BNDRY(3).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (IX.eq.NX) then
      COEF = -RA*P12 + 0.5*RB*(RIJ+R1)
      !      IF (BNDRY(3).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI+HX,YJ)
      IF (BNDRY(4).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (JY.eq.NY) then
      COEF = -Q12 + 0.5*HY*(SIJ+S1)
      !      IF (BNDRY(4).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ+HY)
      if (bndry(5).eq.1) vals(jd) = vals(jd) + coef
    end if
    IF (KZ.eq.NZ) then
      COEF = -RAZ*QB12 + 0.5*RBZ*(SBIJ+SB1)
      !      IF (BNDRY(4).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ+HY)
      if (bndry(6).eq.1) vals(jd) = vals(jd) + coef
    end if

    the_result = 0
    
  end function MATPDE3D_rowFunc


  pure function pc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: pc
    pc = exp(-x*y)
  end function pc

  pure function qc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: qc
    qc = exp(x*y)
  end function qc

  pure function qbc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: qbc
    qbc = exp(x*z)
  end function qbc

  pure function tc(x,y,z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: tc
    tc = 1. / (1.+x+y)
  end function tc

  pure function u(x,y, z)
    real(kind=8), intent(in) :: x, y, z
    real(kind=8) :: u
    u = x * exp(x*y) * sin(pi*x) * sin(pi*y)
  end function u

  pure function rc(beta, x,y, z)
    real(kind=8), intent(in) :: beta, x, y, z
    real(kind=8) :: rc
    rc = beta * (x+y)
  end function rc

  pure function sc(gamma, x,y, z)
    real(kind=8), intent(in) :: gamma, x, y, z
    real(kind=8) :: sc
    sc = gamma * (x+y)
  end function sc

  pure function sbc(gamma, x,y,z)
    real(kind=8), intent(in) :: gamma, x, y, z
    real(kind=8) :: sbc
    sbc = gamma * (y+z)
  end function sbc

  pure function f(beta, gamma, x,y,z)
    real(kind=8), intent(in) :: beta, gamma, x, y, z
    real(kind=8) :: f

    real(kind=8) pxy, pxxy, qxy, qyxy, cxy, rxy, rxxy, sxy, syxy, exy, sx, cx, sy, cy, a, ax, axx, ay, ayy

    pxy  = pc(x,y,z)
    pxxy = -y * exp(-x*y)
    qxy  = qc(x,y,z)
    qyxy = x * exp(x*y)
    cxy  = 1. / ( 1. + x + y)
    rxy  = rc(beta,x,y,z)
    rxxy = beta
    sxy  = sc(gamma,x,y,z)
    syxy = gamma

    exy = exp(x*y)
    sx  = sin(pi*x)
    cx  = cos(pi*x)
    sy  = sin(pi*y)
    cy  = cos(pi*y)

    a   = x * exy * sx * sy
    ax  = exy * (pi * x * cx + (1. + x * y) * sx) * sy
    axx = exy * (pi * (-pi * x * sx + (x * y + 2.) * cx) + y * sx) * sy + y * ax
    ay  = exy * (pi * cy + x * sy) * x * sx
    ayy = x * (exy * pi * (-pi * sy + x * cy) * sx + ay)

    f = - (pxy*axx + qxy*ayy)  + (2.*rxy - pxxy) * ax  +  (2.*sxy - qyxy) * ay  + (rxxy + syxy + cxy) * a

  end function f

end module matpde3d_module
