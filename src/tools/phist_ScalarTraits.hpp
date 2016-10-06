#ifndef PHIST_SCALAR_TRAITS_HPP
#define PHIST_SCALAR_TRAITS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <limits>
#include <complex>
#include <cstdlib>
#include <cstdio>

#include "phist_void_aliases.h"
#include "phist_macros.h"
#include "phist_operator.h"
#include "phist_lapack.h"

#ifdef PHIST_HAVE_GHOST
#include <ghost/types.h>
#endif

#ifdef PHIST_KERNEL_LIB_BUILTIN
extern "C" void phist_Drandom_number(int,double*);
#endif

namespace phist {

template<typename ST>
class ScalarTraits
{
  public: 
  
  typedef typename ST::MissingImplementationOfScalarTraitsClass error;
  
};

#ifdef PHIST_HAVE_SP
template<>
class ScalarTraits< float >
{
  public:
  
  //! alternative typename for ST
  typedef float scalar_t;
#ifdef PHIST_HAVE_GHOST
  static const ghost_datatype ghost_dt = (ghost_datatype)((int)GHOST_DT_FLOAT|(int)GHOST_DT_REAL);
  static const ghost_datatype c_ghost_dt = (ghost_datatype)((int)GHOST_DT_FLOAT|(int)GHOST_DT_COMPLEX);
#endif  
  typedef phist_SlinearOp linearOp_t;
  typedef phist_Smvec mvec_t;
  typedef phist_SsdMat sdMat_t;
  typedef phist_SsparseMat sparseMat_t;
  typedef phist_Sblas_cmplx blas_cmplx_t;
  typedef float blas_scalar_t;
  
  typedef phist_ClinearOp c_linearOp_t; // this is just to allow a simpler implementation of the complex traits class
  typedef phist_Cmvec c_mvec_t; 
  typedef phist_CsparseMat c_sparseMat_t; 
  typedef phist_CsdMat c_sdMat_t; 
  
  //! for complex types, data type of real and imag part.
  //! for real types, magn_t=scalar_t
  typedef float magn_t;

#ifdef PHIST_HAVE_MPI
  //! data type for MPI communication
  static inline MPI_Datatype mpi_type() {return MPI_FLOAT;}
  static inline MPI_Datatype mpi_complex_type() {return MPI_COMPLEX;}
#endif

  //! returns the type prefix as a char, for instance 'S', 'D'
  static inline char type_char()
  {
    return 'S';
  }

  //! returns the type prefix for the corresponding complex data type
  static inline char complex_type_char()
  {
    return 'C';
  }
  
  //! wether ST is a complex data type
  static bool is_complex()
  {
    return false;
  }
    
  //! random number
  static inline scalar_t rand()
  {
    return 2.0f*(float)std::rand()/(float)RAND_MAX-1.0;
  }

  //! same random number on all processed
  static inline scalar_t prand()
  {
    scalar_t tmp = rand();
#ifdef PHIST_HAVE_MPI
    MPI_Bcast(&tmp, 1, mpi_type(), 0, MPI_COMM_WORLD);
#endif
    return tmp;
  }

  //! square root
  static inline scalar_t sqrt(const scalar_t& x)
  {
    return std::sqrt(x);
  }

  //! square root with increased precision. also returns 1/sqrt(x).
  static inline void prec_sqrt(const scalar_t& x, const scalar_t& xC,
                   scalar_t& sqrt, scalar_t& sqrtC,
                   scalar_t& inv_sqrt, scalar_t& inv_sqrtC)
  {
    throw "not implemented";
  }
  
  //! in-place computation of 1/x with increased precision
  static inline void prec_div(const scalar_t& x, const scalar_t& xC,
                               scalar_t& inv,     scalar_t& invC)
  {
    throw "not implemented";
  }
    
  //! absolute value
  static inline magn_t abs(const scalar_t& x)
  {
    return std::abs(x);
  }  
    
  //! complex conjugate
  static inline scalar_t conj(const scalar_t& x)
  {
    return x;
  }
    
  //! real part
  static inline magn_t real(const scalar_t& x)
  {
    return x;
  }
    
  //! imaginary part (0 for real data types)
  static inline magn_t imag(const scalar_t& x)
  {
    PHIST_TOUCH(x)
    return 0.0f;
  }

  //! maximum value
  //! \warning in contrast std::max and fmax actually return NaN when any argument is NaN
  static inline scalar_t max(const scalar_t& x, const scalar_t& y)
  {
    if( x != x )
      return x;
    if( x >  y )
      return x;
    return y;
  }
    
  //! minimum value
  //! \warning in contrast std::max and fmax actually return NaN when any argument is NaN
  static inline scalar_t min(const scalar_t& x, const scalar_t& y)
  {
    if( x != x )
      return x;
    if( x <  y )
      return x;
    return y;
  }
    
  //! scalar 0
  static inline scalar_t zero()
  {
    return 0.0f;
  }
    
  //! scalar 1
  static inline scalar_t one()
  {
    return 1.0f;
  }
    
  //! imaginary unit
  static inline scalar_t cmplx_I()
  {
    return 0.0f;
  }
  
  //! machine epsilon around 1.0 (distance between
  //! 1.0 and the next floating point number)
  static inline magn_t eps(){return std::numeric_limits<magn_t>::epsilon();}

  //! suggested criterion for determining the rank of a matrix.
  
  //! To get uniform behavior we *suggest* to use this value, but it is not compulsory. If s_max is the largest
  //! singular value of a matrix M, a singular value s(M) is considered to be 0 if s_max<rankTol or s<rankTol*s_max.
  //! If the high_prec argument is true (default false), it is asusmed that M is stored in double double precision 
  //! and a value <eps will be returned.
  static inline magn_t rankTol(bool high_prec=false) {return high_prec? 1.0f-12: 10.0f*eps();}

};
#endif
template<>
class ScalarTraits< double >
  {
  public: 
  
  //! alternative typename for ST
  typedef double scalar_t;
#ifdef PHIST_HAVE_GHOST
  static const ghost_datatype ghost_dt = (ghost_datatype)((int)GHOST_DT_DOUBLE|(int)GHOST_DT_REAL);
  static const ghost_datatype c_ghost_dt = (ghost_datatype)((int)GHOST_DT_DOUBLE|(int)GHOST_DT_COMPLEX);
#endif  
  typedef phist_DlinearOp linearOp_t; 
  typedef phist_Dmvec mvec_t; 
  typedef phist_DsparseMat sparseMat_t; 
  typedef phist_DsdMat sdMat_t; 
  typedef phist_Dblas_cmplx blas_cmplx_t;
  typedef double blas_scalar_t;

  typedef phist_ZlinearOp c_linearOp_t; // this is just to allow a simpler implementation of the complex traits class
  typedef phist_Zmvec c_mvec_t;
  typedef phist_ZsparseMat c_sparseMat_t;
  typedef phist_ZsdMat c_sdMat_t;
  
  //! for complex types, data type of real and imag part.
  //! for real types, magn_t=scalar_t
  typedef double magn_t;

#ifdef PHIST_HAVE_MPI
  //! data type for MPI communication
  static inline MPI_Datatype mpi_type() {return MPI_DOUBLE;}
  static inline MPI_Datatype mpi_complex_type() {return MPI_DOUBLE_COMPLEX;}
#endif

  //! returns the type prefix as a char, for instance 'S', 'D'
  static inline char type_char()
  {
    return 'D';
  }

  //! returns the type prefix for the corresponding complex data type
  static inline char complex_type_char()
  {
    return 'Z';
  }
  
  //! wether ST is a complex data type
  static bool is_complex()
  {
    return false;
  }

  //! random number
  static inline scalar_t rand()
  {
    return 2.0*(double)std::rand()/(double)RAND_MAX-1.0;
  }

  //! same random number on all processed
  static inline scalar_t prand()
  {
#ifdef PHIST_KERNEL_LIB_BUILTIN
    double tmp = 0;
    phist_Drandom_number(1,&tmp);
#else
    scalar_t tmp = rand();
#ifdef PHIST_HAVE_MPI
    MPI_Bcast(&tmp, 1, mpi_type(), 0, MPI_COMM_WORLD);
#endif
#endif
    return tmp;
  }

  //! square root
  static inline scalar_t sqrt(const scalar_t& x)
  {
    return std::sqrt(x);
  }

  //! square root with increased precision. also returns 1/sqrt(x).
  static inline void prec_sqrt(const scalar_t& x, const scalar_t& xC,
                   scalar_t& sqrt, scalar_t& sqrtC,
                   scalar_t& inv_sqrt, scalar_t& inv_sqrtC)
  {
    throw "not implemented";
  }
  
  //! in-place computation of 1/x with increased precision
  static inline void prec_div(const scalar_t& x, const scalar_t& xC,
                               scalar_t& inv,     scalar_t& invC)
  {
    throw "not implemented";
  }

  //! absolute value
  static inline magn_t abs(const scalar_t& x)
  {
    return std::abs(x);
  }  
    
  //! complex conjugate
  static inline scalar_t conj(const scalar_t& x)
  {
    return x;
  }
    
  //! real part
  static inline magn_t real(const scalar_t& x)
  {
    return x;
  }
    
  //! imaginary part (0 for real data types)
  static inline magn_t imag(const scalar_t& x)
  {
    PHIST_TOUCH(x)
    return 0.0;
  }

  //! maximum value
  //! \warning in contrast std::max and fmax actually return NaN when any argument is NaN
  static inline scalar_t max(const scalar_t& x, const scalar_t& y)
  {
    if( x != x )
      return x;
    if( x >  y )
      return x;
    return y;
  }
    
  //! minimum value
  //! \warning in contrast std::max and fmax actually return NaN when any argument is NaN
  static inline scalar_t min(const scalar_t& x, const scalar_t& y)
  {
    if( x != x )
      return x;
    if( x <  y )
      return x;
    return y;
  }
    
  //! scalar 0
  static inline scalar_t zero()
  {
    return 0.0;
  }
    
  //! scalar 1
  static inline scalar_t one()
  {
    return 1.0;
  }
    
  //! imaginary unit
  static inline scalar_t cmplx_I()
  {
    return 0.0;
  }
  
  //! machine epsilon around 1.0 (distance between
  //! 1.0 and the next floating point number)
  static inline magn_t eps(){return std::numeric_limits<magn_t>::epsilon();}

  //! suggested criterion for determining the rank of a matrix.
  
  //! To get uniform behavior we *suggest* to use this value, but it is not compulsory. If s_max is the largest
  //! singular value of a matrix M, a singular value s(M) is considered to be 0 if s_max<rankTol or s<rankTol*s_max.
  //! If the high_prec argument is true (default false), it is asusmed that M is stored in double double precision 
  //! and a value <eps will be returned.
  static inline magn_t rankTol(bool high_prec=false) {return high_prec? 1.0e-25: 10.0f*eps();}
};


template<typename MT>
class ScalarTraits< std::complex<MT> >
{
  public: 
  
  //! alternative typename for ST
  typedef typename std::complex<MT> scalar_t;
#ifdef PHIST_HAVE_GHOST
  static const ghost_datatype ghost_dt = ScalarTraits<MT>::c_ghost_dt;
#endif
  typedef typename ScalarTraits<MT>::c_linearOp_t linearOp_t; 
  typedef typename ScalarTraits<MT>::c_mvec_t mvec_t; 
  typedef typename ScalarTraits<MT>::c_sparseMat_t sparseMat_t; 
  typedef typename ScalarTraits<MT>::c_sdMat_t sdMat_t; 
  typedef typename ScalarTraits<MT>::blas_cmplx_t blas_cmplx_t;
  typedef typename ScalarTraits<MT>::blas_cmplx_t blas_scalar_t;

  //! for complex types, data type of real and imag part.
  //! for real types, magn_t=scalar_t
  typedef MT magn_t;

#ifdef PHIST_HAVE_MPI
  
  //! data type for MPI communication. CAVEAT: the MPI data type
  //! is simply an array of two salars, only the MPI_SUM op will
  //! work correctly. Others, like MPI_PROD, will perform the wrong
  //! operation 'result = re(X)*re(Y) + I im(X)*im(Y)'
  static inline MPI_Datatype mpi_type() {return ScalarTraits<magn_t>::mpi_complex_type();}
    
#endif

  //! returns the type prefix as a char, for instance 'S', 'D'
  static inline char type_char()
  {
    return ScalarTraits<MT>::complex_type_char();
  }

  //! wether ST is a complex data type
  static bool is_complex()
  {
    return true;
  }
    
  //! random number
  static inline scalar_t rand()
  {
    return ScalarTraits<magn_t>::rand() + cmplx_I()*ScalarTraits<magn_t>::rand();
  }

  //! same random number on all processed
  static inline scalar_t prand()
  {
    scalar_t tmp = rand();
#ifdef PHIST_HAVE_MPI
    MPI_Bcast(&tmp, 1, mpi_type(), 0, MPI_COMM_WORLD);
#endif
    return tmp;
  }

  //! square root
  static inline scalar_t sqrt(const scalar_t& x)
  {
    return std::sqrt(x);
  }

  //! square root with increased precision. also returns 1/sqrt(x).
  static inline void prec_sqrt(const scalar_t& x, const scalar_t& xC,
                   scalar_t& sqrt, scalar_t& sqrtC,
                   scalar_t& inv_sqrt, scalar_t& inv_sqrtC)
  {
    throw "not implemented";
  }
  
  //! in-place computation of 1/x with increased precision
  static inline void prec_div(const scalar_t& x, const scalar_t& xC,
                               scalar_t& inv,     scalar_t& invC)
  {
    throw "not implemented";
  }
    
  //! absolute value
  static inline magn_t abs(const scalar_t& x)
  {
    return std::abs(x);
  }  
    
  //! complex conjugate
  static inline scalar_t conj(const scalar_t& x)
  {
    return std::conj(x);
  }
    
  //! real part
  static inline magn_t real(const scalar_t& x)
  {
    return std::real(x);
  }
    
  //! imaginary part (0 for real data types)
  static inline magn_t imag(const scalar_t& x)
  {
    return std::imag(x);
  }
    
  //! maximum value
  //! \warning in contrast std::max and fmax actually return NaN when any argument is NaN
  static inline scalar_t max(const scalar_t& x, const scalar_t& y)
  {
    if( x != x )
      return x;
    if( std::abs(x) >  std::abs(y) )
      return x;
    return y;
  }
    
  //! minimum value
  //! \warning in contrast std::max and fmax actually return NaN when any argument is NaN
  static inline scalar_t min(const scalar_t& x, const scalar_t& y)
  {
    if( x != x )
      return x;
    if( std::abs(x) <  std::abs(y) )
      return x;
    return y;
  }
    
  //! scalar 0
  static inline scalar_t zero()
  {
    return scalar_t(ScalarTraits<MT>::zero(),ScalarTraits<MT>::zero());
  }
    
  //! scalar 1
  static inline scalar_t one()
  {
    return scalar_t(ScalarTraits<MT>::one(),ScalarTraits<MT>::zero());
  }
    
  //! imaginary unit
  static inline scalar_t cmplx_I()
  {
    return scalar_t(ScalarTraits<MT>::zero(),ScalarTraits<MT>::one());
  }
  
  //! machine epsilon around 1.0 (distance between
  //! 1.0 and the next floating point number)
  static inline magn_t eps(){return ScalarTraits<MT>::eps();}

  //! suggested criterion for determining the rank of a matrix.
  
  //! To get uniform behavior we *suggest* to use this value, but it is not compulsory. If s_max is the largest
  //! singular value of a matrix M, a singular value s(M) is considered to be 0 if s_max<rankTol or s<rankTol*s_max.
  //! If the high_prec argument is true (default false), it is asusmed that M is stored in double double precision 
  //! and a value <eps will be returned.
  static inline magn_t rankTol(bool high_prec=false) {return ScalarTraits<MT>::rankTol(high_prec);}
};

  
}//namepsace



#endif
