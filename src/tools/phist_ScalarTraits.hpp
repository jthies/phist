#ifndef PHIST_SCALAR_TRAITS_HPP
#define PHIST_SCALAR_TRAITS_HPP

#include <limits>
#include <complex>
#include <cstdlib>

#include "phist_void_aliases.h"
#include "phist_operator.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_HAVE_GHOST
#include "ghost.h"
#endif

namespace phist {

template<typename ST>
class ScalarTraits
  {
  public: 
  
  // this should prevent compiling a non-specialized
  // template instantiation, but I'm not sure if it does.
  typedef typename ST::missing_specialization_for_class_ScalarTraits bar;
  };

template<>
class ScalarTraits< float >
  {
  public:
  
  //! alternative typename for ST
  typedef float scalar_t;
#ifdef PHIST_HAVE_GHOST
  static const int ghost_dt = GHOST_BINCRS_DT_FLOAT|GHOST_BINCRS_DT_REAL;
  static const int c_ghost_dt = GHOST_BINCRS_DT_FLOAT|GHOST_BINCRS_DT_COMPLEX;
#endif  
  typedef Sop_t op_t;
  typedef Smvec_t mvec_t;
  typedef SsdMat_t sdMat_t;
  typedef ScrsMat_t crsMat_t;
  
  typedef Cop_t c_op_t; // this is just to allow a simpler implementation of the complex traits class
  typedef Cmvec_t c_mvec_t; 
  typedef CcrsMat_t c_crsMat_t; 
  typedef CsdMat_t c_sdMat_t; 
  
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
  static inline const char complex_type_char()
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

  //! square root
  static inline scalar_t sqrt(const scalar_t& x)
    {
    return std::sqrt(x);
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
    return 0.0f;
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
  static inline scalar_t I()
    {
    return 0.0f;
    }
  
  //! machine epsilon around 1.0 (distance between
  //! 1.0 and the next floating point number)
  static inline magn_t eps(){return std::numeric_limits<magn_t>::epsilon();}
  };

template<>
class ScalarTraits< double >
  {
  public: 
  
  //! alternative typename for ST
  typedef double scalar_t;
#ifdef PHIST_HAVE_GHOST
  static const int ghost_dt = GHOST_BINCRS_DT_DOUBLE|GHOST_BINCRS_DT_REAL;
  static const int c_ghost_dt = GHOST_BINCRS_DT_DOUBLE|GHOST_BINCRS_DT_COMPLEX;
#endif  
  typedef Dop_t op_t; 
  typedef Dmvec_t mvec_t; 
  typedef DcrsMat_t crsMat_t; 
  typedef DsdMat_t sdMat_t; 
  
  typedef Zop_t c_op_t; // this is just to allow a simpler implementation of the complex traits class
  typedef Zmvec_t c_mvec_t;
  typedef ZcrsMat_t c_crsMat_t;
  typedef ZsdMat_t c_sdMat_t;
  
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
  static inline const char complex_type_char()
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
    
  //! square root
  static inline scalar_t sqrt(const scalar_t& x)
    {
    return std::sqrt(x);
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
    return 0.0;
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
  static inline scalar_t I()
    {
    return 0.0;
    }
  
  //! machine epsilon around 1.0 (distance between
  //! 1.0 and the next floating point number)
  static inline magn_t eps(){return std::numeric_limits<magn_t>::epsilon();}
  };


template<typename MT>
class ScalarTraits< std::complex<MT> >
  {
  public: 
  
  //! alternative typename for ST
  typedef typename std::complex<MT> scalar_t;
#ifdef PHIST_HAVE_GHOST
  static const int ghost_dt = ScalarTraits<MT>::c_ghost_dt;
#endif
  typedef typename ScalarTraits<MT>::c_op_t op_t; 
  typedef typename ScalarTraits<MT>::c_mvec_t mvec_t; 
  typedef typename ScalarTraits<MT>::c_crsMat_t crsMat_t; 
  typedef typename ScalarTraits<MT>::c_sdMat_t sdMat_t; 

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
    return ScalarTraits<magn_t>::rand() + I()*ScalarTraits<magn_t>::rand();
    }

  //! square root
  static inline scalar_t sqrt(const scalar_t& x)
    {
    return std::sqrt(x);
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
  static inline scalar_t I()
    {
    return scalar_t(ScalarTraits<MT>::zero(),ScalarTraits<MT>::one());
    }
  
  //! machine epsilon around 1.0 (distance between
  //! 1.0 and the next floating point number)
  static inline magn_t eps(){return ScalarTraits<MT>::eps();}
  };

  
}//namepsace



#endif
