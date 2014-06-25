! gemm kernels for mvec_times_sdmat

subroutine dgemm_sB_1_k(nrows,k,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(1,k)
  real(kind=8), intent(inout) :: C(1,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_1_k

subroutine dgemm_sB_2_k(nrows,k,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(2,k)
  real(kind=8), intent(inout) :: C(2,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_2_k

subroutine dgemm_sB_4_k(nrows,k,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(4,k)
  real(kind=8), intent(inout) :: C(4,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_4_k

subroutine dgemm_sB_8_k(nrows,k,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(k,nrows)
  real(kind=8), intent(in) :: B(8,k)
  real(kind=8), intent(inout) :: C(8,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_8_k



subroutine dgemm_sB_k_1(nrows,n,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, n
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(1,nrows)
  real(kind=8), intent(in) :: B(1,n)
  real(kind=8), intent(inout) :: C(n,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_k_1

subroutine dgemm_sB_k_2(nrows,n,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, n
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(2,nrows)
  real(kind=8), intent(in) :: B(2,n)
  real(kind=8), intent(inout) :: C(n,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_k_2

subroutine dgemm_sB_k_4(nrows,n,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, n
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(4,nrows)
  real(kind=8), intent(in) :: B(4,n)
  real(kind=8), intent(inout) :: C(n,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_k_4


subroutine dgemm_sB_k_8(nrows,n,alpha, A, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, n
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(8,nrows)
  real(kind=8), intent(in) :: B(8,n)
  real(kind=8), intent(inout) :: C(n,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_k_8


subroutine dgemm_sB_strided_1_k(nrows,k,alpha, A, lda, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(1,k)
  real(kind=8), intent(inout) :: C(1,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_strided_1_k

subroutine dgemm_sB_strided_2_k(nrows,k,alpha, A, lda, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(2,k)
  real(kind=8), intent(inout) :: C(2,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_strided_2_k

subroutine dgemm_sB_strided_4_k(nrows,k,alpha, A, lda, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(4,k)
  real(kind=8), intent(inout) :: C(4,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_strided_4_k

subroutine dgemm_sB_strided_8_k(nrows,k,alpha, A, lda, B, beta, C)
  implicit none
  integer, intent(in) :: nrows, k, lda
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(8,k)
  real(kind=8), intent(inout) :: C(8,nrows)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = 0.
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      C(:,i) = beta*C(:,i)
      do j = 1, k, 1
        C(:,i) = C(:,i) + alpha*B(:,j)*A(j,i)
      end do
    end do

  end if
end subroutine dgemm_sB_strided_8_k



subroutine dgemm_sB_strided_k_1(nrows,n,alpha, A, B, beta, C, ldC)
  implicit none
  integer, intent(in) :: nrows, n, ldc
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(1,nrows)
  real(kind=8), intent(in) :: B(1,n)
  real(kind=8), intent(inout) :: C(ldC,*)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_strided_k_1

subroutine dgemm_sB_strided_k_2(nrows,n,alpha, A, B, beta, C, ldC)
  implicit none
  integer, intent(in) :: nrows, n, ldc
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(2,nrows)
  real(kind=8), intent(in) :: B(2,n)
  real(kind=8), intent(inout) :: C(ldC,*)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_strided_k_2

subroutine dgemm_sB_strided_k_4(nrows,n,alpha, A, B, beta, C, ldC)
  implicit none
  integer, intent(in) :: nrows, n, ldc
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(4,nrows)
  real(kind=8), intent(in) :: B(4,n)
  real(kind=8), intent(inout) :: C(ldC,*)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_strided_k_4


subroutine dgemm_sB_strided_k_8(nrows,n,alpha, A, B, beta, C, ldC)
  implicit none
  integer, intent(in) :: nrows, n, ldc
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(8,nrows)
  real(kind=8), intent(in) :: B(8,n)
  real(kind=8), intent(inout) :: C(ldC,*)
  integer :: i, j

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(:,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, nrows, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(:,i))
      end do
    end do

  end if
end subroutine dgemm_sB_strided_k_8


subroutine dgemm_sB_generic(m,n,k,alpha,A,lda,B,beta,C,ldc)
  implicit none
  integer, intent(in) :: m,n,k,lda,ldc
  real(kind=8), intent(in) :: alpha, beta
  real(kind=8), intent(in) :: A(lda,*)
  real(kind=8), intent(in) :: B(k,n)
  real(kind=8), intent(inout) :: C(ldc,*)
  integer :: i, j


  !call dgemm('T','N',n,m,k,alpha,B,k,A,lda,beta,C,ldc)

  if( beta .eq. 0 ) then
!$omp parallel do schedule(static)
    do i = 1, m, 1
      do j = 1, n, 1
        C(j,i) = alpha*sum(B(:,j)*A(1:k,i))
      end do
    end do
  else
!$omp parallel do schedule(static)
    do i = 1, m, 1
      do j = 1, n, 1
        C(j,i) = beta*C(j,i) + alpha*sum(B(:,j)*A(1:k,i))
      end do
    end do
  end if

end subroutine dgemm_sB_generic


subroutine dgemm_sB_generic_inplace(m,n,k,A,lda,B)
  implicit none
  integer, intent(in) :: m,n,k,lda
  real(kind=8), intent(inout) :: A(lda,*)
  real(kind=8), intent(in) :: B(k,n)
  integer :: i, j
  real(kind=8) :: work(n)

!$omp parallel do private(work) schedule(static)
  do i = 1, m, 1
    do j = 1, n, 1
      work(j) = sum(A(1:k,i)*B(:,j))
    end do
    A(1:n,i) = work(:)
  end do

end subroutine dgemm_sB_generic_inplace

