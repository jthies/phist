# try to detect if we can use certain C++11 and Fortran/OpenMP features
# sets the following variables:
# PHIST_HAVE_CXX11_LAMBDAS
# PHIST_HAVE_CXX11_MOVEDEFAULT
# PHIST_HAVE_CXX11_THREADLOCAL
# PHIST_HAVE_OPENMP_SIMD

include (CMakePushCheckState)
include(CheckCXXSourceCompiles)
include(CheckFortranSourceCompiles)

cmake_push_check_state() # Save variables

set(CMAKE_REQUIRED_FLAGS "-std=c++11")

CHECK_CXX_SOURCE_COMPILES("
  int main(int argc, char** argv)
  {
    auto lfunc = [&](){ int _argc = argc; };
    return 0;
  }
" PHIST_HAVE_CXX11_LAMBDAS
)
CHECK_CXX_SOURCE_COMPILES("
  int main(int argc, char** argv)
  {
    thread_local int x;
    return 0;
  }
" PHIST_HAVE_CXX11_THREADLOCAL)

CHECK_CXX_SOURCE_COMPILES("
  template<typename T>
  class C {
  public:
  C(){}
  C(C&&)=default;
  };
  int main(int argc, char** argv)
  {
    C<int> c;
  }
" PHIST_HAVE_CXX11_MOVEDEFAULT)

CHECK_Fortran_SOURCE_COMPILES("
  program test
    implicit none
    integer :: i
    real :: x(8)
    !$omp simd
    do i = 1, 8
      x(i) = 0.
    end do
  end program test
" PHIST_HAVE_OPENMP_SIMD
)

cmake_pop_check_state() # Recover variables
