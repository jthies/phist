!> \file gemm_kernels_sB_augmented.f90
!! Fast parallel BLAS-gemm like subroutines augmented with the product of the result for different blocksizes for mvec_times_sdmat_augmented in mvec_module
!! \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
!!

subroutine dgemm_sB_augmented_1_k(nrows,k,alpha, A, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(1,k)
  real(kind=8), intent(inout) :: C(1,nrows)
  real(kind=8), intent(out) :: D(1,1)
  integer :: i, j
  real(kind=8) :: Dtmp
!dir$ assume_aligned A:64, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      Dtmp = Dtmp + C(1,i)*C(1,i)
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      Dtmp = Dtmp + C(1,i)*C(1,i)
    end do

  end if

  D = Dtmp
end subroutine dgemm_sB_augmented_1_k

subroutine dgemm_sB_augmented_2_k(nrows,k,alpha, A, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(2,k)
  real(kind=8), intent(inout) :: C(2,nrows)
  real(kind=8), intent(out) :: D(2,2)
  integer :: i, j
  real(kind=8) :: Dtmp(2,2)
!dir$ assume_aligned A:64, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 2, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 2, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do

  end if

  ! compose result
  do j = 1, 2, 1
    do i = j, 2, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do
end subroutine dgemm_sB_augmented_2_k

subroutine dgemm_sB_augmented_4_k(nrows,k,alpha, A, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(4,k)
  real(kind=8), intent(inout) :: C(4,nrows)
  real(kind=8), intent(out) :: D(4,4)
  integer :: i, j
  real(kind=8) :: Dtmp(4,4)
!dir$ assume_aligned A:64, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 4, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 4, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do

  end if

  ! compose result
  do j = 1, 4, 1
    do i = j, 4, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do
end subroutine dgemm_sB_augmented_4_k

subroutine dgemm_sB_augmented_8_k(nrows,k,alpha, A, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(8,k)
  real(kind=8), intent(inout) :: C(8,nrows)
  real(kind=8), intent(out) :: D(8,8)
  integer :: i, j
  real(kind=8) :: Dtmp(8,8)
!dir$ assume_aligned A:64, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 8, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 8, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do

  end if

  ! compose result
  do j = 1, 8, 1
    do i = j, 8, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do
end subroutine dgemm_sB_augmented_8_k



subroutine dgemm_sB_augmented_1_strided_k(nrows,k,alpha, A, lda, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(1,k)
  real(kind=8), intent(inout) :: C(1,nrows)
  real(kind=8), intent(out) :: D(1,1)
  integer :: i, j
  real(kind=8) :: Dtmp
!dir$ assume_aligned A:8, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      Dtmp = Dtmp + C(1,i)*C(1,i)
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      Dtmp = Dtmp + C(1,i)*C(1,i)
    end do

  end if

  D = Dtmp
end subroutine dgemm_sB_augmented_1_strided_k

subroutine dgemm_sB_augmented_2_strided_k(nrows,k,alpha, A, lda, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(2,k)
  real(kind=8), intent(inout) :: C(2,nrows)
  real(kind=8), intent(out) :: D(2,2)
  integer :: i, j
  real(kind=8) :: Dtmp(2,2)
!dir$ assume_aligned A:8, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 2, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 2, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do

  end if

  ! compose result
  do j = 1, 2, 1
    do i = j, 2, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do
end subroutine dgemm_sB_augmented_2_strided_k

subroutine dgemm_sB_augmented_4_strided_k(nrows,k,alpha, A, lda, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(4,k)
  real(kind=8), intent(inout) :: C(4,nrows)
  real(kind=8), intent(out) :: D(4,4)
  integer :: i, j
  real(kind=8) :: Dtmp(4,4)
!dir$ assume_aligned A:8, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 4, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 4, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do

  end if

  ! compose result
  do j = 1, 4, 1
    do i = j, 4, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do
end subroutine dgemm_sB_augmented_4_strided_k

subroutine dgemm_sB_augmented_8_strided_k(nrows,k,alpha, A, lda, B, beta, C, D)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(8,k)
  real(kind=8), intent(inout) :: C(8,nrows)
  real(kind=8), intent(out) :: D(8,8)
  integer :: i, j
  real(kind=8) :: Dtmp(8,8)
!dir$ assume_aligned A:8, B:64, C:64

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
      do j = 1, 8, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if

  ! compose result
  do j = 1, 8, 1
    do i = j, 8, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do

end subroutine dgemm_sB_augmented_8_strided_k



subroutine dgemm_sB_augmented_generic(m,n,k,alpha,A,lda,B,beta,C,ldc, D)
  implicit none
  integer, intent(in) :: m,n,k,lda,ldc
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(k,n)
  real(kind=8), intent(inout) :: C(ldc,*)
  real(kind=8), intent(out) :: D(n,n)
  integer :: i, j
  real(kind=8) :: Dtmp(n,n)
!dir$ assume_aligned A:8, B:64, C:8


  !call dgemm('T','N',n,m,k,alpha,B,k,A,lda,beta,C,ldc)

  if( beta .eq. 0 ) then
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, m, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(1:k,i))
      end do
      do j = 1, n, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:n,i)
      end do
    end do
  else
    Dtmp = 0.
!$omp parallel do reduction(+:Dtmp) schedule(static)
    do i = 1, m, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(1:k,i))
      end do
      do j = 1, n, 1
        Dtmp(j:,j) = Dtmp(j:,j) + C(j,i)*C(j:n,i)
      end do
    end do
  end if

  ! compose result
  do j = 1, n, 1
    do i = j, n, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do

end subroutine dgemm_sB_augmented_generic


subroutine dgemm_sB_augmented_generic_inplace(m,n,k,A,lda,B, D)
  implicit none
  integer, intent(in) :: m,n,k,lda
  real(kind=8), intent(inout) :: A(lda,*)
  real(kind=8), intent(in) :: B(k,n)
  real(kind=8), intent(out) :: D(n,n)
  integer :: i, j
  real(kind=8) :: Dtmp(n,n)
  real(kind=8) :: work(n)
!dir$ assume_aligned A:8, B:64

    Dtmp = 0.
!$omp parallel do private(work) reduction(+:Dtmp) schedule(static)
  do i = 1, m, 1
    do j = 1, n, 1
      work(j) = sum(A(1:k,i)*B(:,j))
    end do
    do j = 1, n, 1
      Dtmp(j:,j) = Dtmp(j:,j) + work(j)*work(j:)
    end do
    A(1:n,i) = work(:)
  end do

  ! compose result
  do j = 1, n, 1
    do i = j, n, 1
      D(i,j) = Dtmp(i,j)
      D(j,i) = Dtmp(i,j)
    end do
  end do

end subroutine dgemm_sB_augmented_generic_inplace

