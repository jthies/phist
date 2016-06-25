#include "phist_config.h"
#ifdef PHIST_HAVE_GHOST
#include "ghost/config.h"
#endif

#ifndef PHIST_HAVE_GHOST
#define G_LIDX_T C_INT32_T
# ifndef PHIST_FORCE_INT_GIDX
# define G_GIDX_T C_INT32_T
# else
# define G_GIDX_T C_INT64_T
# endif
#else
#ifdef GHOST_IDX64_LOCAL
#define G_LIDX_T C_INT64_T
#else
#define G_LIDX_T C_INT32_T
#endif
#ifdef GHOST_IDX64_GLOBAL
#define G_GIDX_T C_INT64_T
#else
#define G_GIDX_T C_INT32_T
#endif
#endif

! From the NEP matrix collection in the matrix market MATPDE matrix generator
module matpde_module
  implicit none

  public :: MATPDE_rowFunc
  public :: MATPDE_initDimensions

  real(kind=8), parameter :: pi = 4*atan(1.0_8)

  integer(kind=8) :: nx, ny
  integer :: level

contains

  subroutine MATPDE_initDimensions(new_nx, new_ny, nrows, maxnne_per_row) bind(C, name='MATPDE_initDimensions')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: new_nx, new_ny
    integer(kind=G_GIDX_T), intent(out) :: nrows
    integer(kind=G_LIDX_T), intent(out) :: maxnne_per_row

    if( new_nx .ne. new_ny ) then
      write(*,*) 'MATPDE: error, nx != ny!'
      call exit(1)
    end if
    level = nint(log(1._8*new_nx)/log(2._8))
    if( 2**level .ne. new_nx ) then
      write(*,*) 'MATPDE: error, nx needs to be a power of 2!'
      call exit(1)
    end if

    nx = new_nx
    ny = new_ny

    nrows = nx*ny

    maxnne_per_row = 5

  end subroutine MATPDE_initDimensions


  ! useful quadtree ordering
  pure function idOfCoord(coord) result(id)
    integer(kind=8), intent(in) :: coord(2)
    integer(kind=8) :: id
    integer(kind=8) :: fak2, fak4
    integer(kind=8) :: upperBound
    integer(kind=8) :: myCoord(2)
    integer :: i

    upperBound = 2**level
    myCoord = modulo(coord, upperBound)

    fak2 = 1
    fak4 = 1
    id = 1
    do i = 1, level
      id = id + fak4 * ( 1 * mod(myCoord(1)/fak2,2_8) &
        &              + 2 * mod(myCoord(2)/fak2,2_8) )
      fak2 = fak2 * 2
      fak4 = fak4 * 4
    end do
  end function idOfCoord


  pure function coordOfId(id) result(coord)
  
    use, intrinsic :: iso_c_binding
  
    integer(kind=G_GIDX_T), intent(in) :: id
    integer(kind=8) :: coord(2)
    integer(kind=8) :: tElem, fak(2)
    integer :: bitlevel, i

    fak(1) = 1
    fak(2) = 2

    bitlevel = 1
    tElem = mod(id - 1, 4**level)

    coord = 0
    do while( (tElem / fak(1)) > 0 )
      do i = 1, 2
        coord(i) = coord(i) + bitlevel * int(mod(tElem/fak(i),2_8))
      end do
      bitlevel = bitlevel * 2
      fak = fak * 4
    end do
  end function coordOfId


  function MATPDE_rowFunc(row, nnz, cols, vals, dummy) result(the_result) bind(C, name='MATPDE_rowFunc')
    use, intrinsic :: iso_c_binding
    integer(G_GIDX_T), value :: row
    integer(G_LIDX_T), intent(inout) :: nnz
    integer(G_GIDX_T), intent(inout) :: cols(*)
    real(C_DOUBLE),    intent(inout) :: vals(*)
    TYPE(c_ptr), value               :: dummy
    integer(C_INT) :: the_result

! copied from MatrixMarket, original description:
    !
    !  Purpose
    !  =======
    !
    !  Forms the five-point central differences operator for the elliptic
    !  PDE
    !
    !   -(P Ux)x -(Q Uy)y  +  R Ux + (R U)x  +  S Uy + (S U)y + T U = F
    !
    !  where P, Q, R, S and T are the functions of x and y.  The domain
    !  is the unit square (0,1)x(0,1), and the boundary condition are
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
    real(kind=8), parameter :: beta = 20., gamma = 0.
    !     ..
    !     .. Scalar variables ..
    integer(kind=8) :: ix, jy, index, jd, coord(2), coord_(2)
    real(kind=8) :: bndry( 4 )
    real(kind=8) :: hx, hy, hy2, ra, rb, coef
    real(kind=8) :: p12, pm12, q12, qm12, xi, rij, r1, rm1, sij, s1, sm1, yj

    !     set up initial boundary conditions
    BNDRY( 1 ) = 0.
    BNDRY( 2 ) = 0.
    BNDRY( 3 ) = 0.
    BNDRY( 4 ) = 0.

    !     set up other parameters
    HX = 1. / (DBLE(NX) + 1.)
    HY = 1. / (DBLE(NY) + 1.)
    HY2 = HY*HY
    RA = (HY/HX)**2
    RB  = HY2/HX

    nnz = 0
    ! do jy = 1, ny
    !   do ix = 1, nx
    coord = coordOfId(row+1)
    ix = coord(1)+1
    jy = coord(2)+1
    !jy = row/nx+1
    !ix = mod(row,nx)+1
      
    YJ = HY*DBLE(JY)
    XI = HX*DBLE(IX)
    !INDEX = IX + (JY-1)*NX
    !    RHS(INDEX) = HY2*F(BETA,GAMMA,XI,YJ)
    P12  = PC (XI + 0.5*HX,YJ)
    PM12 = PC (XI - 0.5*HX,YJ)
    Q12  = QC (XI,YJ + 0.5*HY)
    QM12 = QC (XI,YJ - 0.5*HY)
    RIJ  = RC (BETA,XI,YJ)
    R1   = RC (BETA,XI + HX,YJ)
    RM1  = RC (BETA,XI - HX,YJ)
    SIJ  = SC (GAMMA,XI,YJ)
    S1   = SC (GAMMA,XI,YJ + HY)
    SM1  = SC (GAMMA,XI,YJ - HY)

    !           DIAGONAL.
    nnz = nnz + 1
    vals(nnz) = ra*(p12 + pm12)  +  q12 + qm12  + hy2*tc(xi,yj)
    cols(nnz) = idOfCoord(coord) !index
    jd = nnz

    !           LOWEST BAND.
    if (jy.ne.1) then
      nnz = nnz + 1
      vals(nnz) = - ( qm12 + 0.5*hy*(sij+sm1) )
      coord_ = coord-(/0,1/)
      cols(nnz) = idOfCoord(coord_) !index - nx
    end if

    !           SUB-DIAGONAL.
    if (ix.ne.1) then
      nnz = nnz + 1
      vals(nnz) = - ( ra*pm12 + 0.5*rb*(rij+rm1) )
      coord_ = coord-(/1,0/)
      cols(nnz) = idOfCoord(coord_) !index - 1
    end if

    !           SUPER-DIAGONAL.
    if (ix.ne.nx)  then
      nnz = nnz + 1
      vals(nnz) = -ra*p12 + 0.5*rb*(rij+r1)
      coord_ = coord+(/1,0/)
      cols(nnz) = idOfCoord(coord_) !index + 1
    end if

    !           HIGHEST BAND.
    if (jy.ne.ny) then
      nnz = nnz + 1
      vals(nnz) = -q12 + 0.5*hy*(sij+s1)
      coord_ = coord+(/0,1/)
      cols(nnz) = idOfCoord(coord_) !index + nx
    end if

    !           BOUNDARY CONDITIONS.
    if (jy.eq.1) then
      COEF = - ( QM12 + 0.5*HY*(SIJ+SM1) )
      !      IF (BNDRY(1).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ-HY)
      IF (BNDRY(1).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    if (ix.eq.1) then
      COEF = - ( RA*PM12 + 0.5*RB*(RIJ+RM1) )
      !      IF (BNDRY(2).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI-HX,YJ)
      IF (BNDRY(2).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (IX.eq.NX) then
      COEF = -RA*P12 + 0.5*RB*(RIJ+R1)
      !      IF (BNDRY(3).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI+HX,YJ)
      IF (BNDRY(3).EQ.1) vals(JD) = vals(JD) + COEF
    end if
    IF (JY.eq.NY) then
      COEF = -Q12 + 0.5*HY*(SIJ+S1)
      !      IF (BNDRY(4).EQ.0) RHS(INDEX) = RHS(INDEX) - COEF*U(XI,YJ+HY)
      if (bndry(4).eq.1) vals(jd) = vals(jd) + coef
    end if

    cols(1:nnz) = cols(1:nnz) - 1

    the_result = 0
    
  end function MATPDE_rowFunc


  pure function pc(x,y)
    real(kind=8), intent(in) :: x, y
    real(kind=8) :: pc
    pc = exp(-x*y)
  end function pc

  pure function qc(x,y)
    real(kind=8), intent(in) :: x, y
    real(kind=8) :: qc
    qc = exp(x*y)
  end function qc

  pure function tc(x,y)
    real(kind=8), intent(in) :: x, y
    real(kind=8) :: tc
    tc = 1. / (1.+x+y)
  end function tc

  pure function u(x,y)
    real(kind=8), intent(in) :: x, y
    real(kind=8) :: u
    u = x * exp(x*y) * sin(pi*x) * sin(pi*y)
  end function u

  pure function rc(beta, x,y)
    real(kind=8), intent(in) :: beta, x, y
    real(kind=8) :: rc
    rc = beta * (x+y)
  end function rc

  pure function sc(gamma, x,y)
    real(kind=8), intent(in) :: gamma, x, y
    real(kind=8) :: sc
    sc = gamma * (x+y)
  end function sc

  pure function f(beta, gamma, x,y)
    real(kind=8), intent(in) :: beta, gamma, x, y
    real(kind=8) :: f

    real(kind=8) pxy, pxxy, qxy, qyxy, cxy, rxy, rxxy, sxy, syxy, exy, sx, cx, sy, cy, a, ax, axx, ay, ayy

    pxy  = pc(x,y)
    pxxy = -y * exp(-x*y)
    qxy  = qc(x,y)
    qyxy = x * exp(x*y)
    cxy  = 1. / ( 1. + x + y)
    rxy  = rc(beta,x,y)
    rxxy = beta
    sxy  = sc(gamma,x,y)
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

end module matpde_module
