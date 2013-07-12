#ifndef PHIST_SCALAR_TRAITS_HPP
#define PHIST_SCALAR_TRAITS_HPP

#include <limits>
#include <complex>

namespace phist {

template<typename ST>
class ScalarTraits
  {
  public: 
  
  // if no specialization exists, this class can't do anything
  // and a compile-time error should occur
  typedef ST scalar_t;
//  typedef missing_specialization_for_class_ScalarTraits scalar_t;
  };

template<>
class ScalarTraits< float >
  {
  public:
  
  //! alternative typename for ST
  typedef float scalar_t;
  //! for complex types, data type of real and imag part.
  //! for real types, magn_t=scalar_t
  typedef float magn_t;

  //! returns the type prefix as a char, for instance 'S', 'D'
  static inline const char type_char()
    {
    return 'S';
    }

  //! returns the type prefix for the corresponding complex data type
  static inline const char c_type_char()
    {
    return 'C';
    }
  
  //! wether ST is a complex data type
  static bool is_complex()
    {
    return false;
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
  //! for complex types, data type of real and imag part.
  //! for real types, magn_t=scalar_t
  typedef double magn_t;

  //! returns the type prefix as a char, for instance 'S', 'D'
  static inline const char type_char()
    {
    return 'D';
    }

  //! returns the type prefix for the corresponding complex data type
  static inline const char c_type_char()
    {
    return 'Z';
    }
  
  //! wether ST is a complex data type
  static bool is_complex()
    {
    return false;
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
  //! for complex types, data type of real and imag part.
  //! for real types, magn_t=scalar_t
  typedef MT magn_t;

  //! returns the type prefix as a char, for instance 'S', 'D'
  static inline const char type_char()
    {
    return ScalarTraits<MT>::complex_type_char();
    }

  //! wether ST is a complex data type
  static bool is_complex()
    {
    return true;
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
