# simple script to test for SSE, AVX, AVX2 and AVX512.
# currently does not support cross-compiling...
#
# on output, PHIST_HAVE_[SSE,AVX,AVX2,AVX512] are defined (or not)

include(CheckCSourceRuns)

set (CMAKE_C_FLAGS        "-march=native")
  
CHECK_C_SOURCE_RUNS("
#include \"emmintrin.h\"
#include \"immintrin.h\"
  int main(){
  double a=0.0;
  __m128d a_ = _mm_set_pd(a,a);\
  return 0;
  }
" PHIST_HAVE_SSE)

CHECK_C_SOURCE_RUNS("
#include \"emmintrin.h\"
#include \"immintrin.h\"
int main(){
double a=0.0;
__m256d a_ = _mm256_set_pd(a,a,a,a);\
return 0;
}
" PHIST_HAVE_AVX)

CHECK_C_SOURCE_RUNS("
#include \"emmintrin.h\"
#include \"immintrin.h\"
int main(){
double a=1.0,b=2.0,c=3.0;
__m256d a_ = _mm256_set_pd(a,a,a,a);\
__m256d b_ = _mm256_set_pd(b,b,b,b);\
__m256d c_ = _mm256_set_pd(c,c,c,c);\
__m256d d_  = _mm256_fmadd_pd(a_,b_,c_);
return 0;
}
" PHIST_HAVE_AVX2)

CHECK_C_SOURCE_RUNS("
#include \"emmintrin.h\"
#include \"immintrin.h\"
int main(){
double a=0.0;
__m512d a_ = _mm512_set_pd(a,a,a,a,a,a,a,a);\
return 0;
}
" PHIST_HAVE_AVX512)
