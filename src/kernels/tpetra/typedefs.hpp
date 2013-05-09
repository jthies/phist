#ifndef KERNELS_TPETRA_TYPEDEFS_HPP
#define KERNELS_TPETRA_TYPEDEFS_HPP

#include "Teuchos_Comm.hpp"

#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"

typedef int lidx_t; // local ordinal/index type
typedef int gidx_t; // global ordinal/index type


typedef Kokkos::DefaultNode::DefaultNodeType node_t; // from the Kokkos node API
typedef Tpetra::Map<lidx_t,gidx_t,node_t> map_t;
typedef Teuchos::Comm<int> comm_t;

template<typename ST>
class Traits
  {
  public:
  
  //!
  typedef typename Kokkos::DefaultKernels<ST,lidx_t,node_t>::SparseOps localOps_t;

  //! multi vectors
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> mvec_t;

  //! serial dense matrix - just a multivector with a serial map.
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> sdMat_t;

  //! CRS matrices
  typedef Tpetra::CrsMatrix<ST,lidx_t,gidx_t,node_t,localOps_t> crsMat_t;

  //! for performing the MVM
  typedef Tpetra::CrsMatrixMultiplyOp<ST,ST,lidx_t,gidx_t,node_t,localOps_t> crsMVM_t;

  //! scalar 1
  static inline ST one(){return Teuchos::ScalarTraits<ST>::one();}

  //! scalar 0
  static inline ST zero(){return Teuchos::ScalarTraits<ST>::zero();}

  
  };



#ifdef _IS_COMPLEX_
namespace Teuchos {
// Partial specialization for std::complex numbers templated on real type <_MT_>
struct ScalarTraits<_ST_>
{
  typedef _ST_  ComplexT;
  typedef std::complex<typename ScalarTraits<_MT_>::halfPrecision> halfPrecision;
  typedef std::complex<typename ScalarTraits<_MT_>::doublePrecision> doublePrecision;
  typedef typename ScalarTraits<<_MT_>>::magnitudeType magnitudeType;
  static const bool isComplex = true;
  static const bool isOrdinal = ScalarTraits<<_MT_>>::isOrdinal;
  static const bool isComparable = false;
  static const bool hasMachineParameters = true;
  static inline magnitudeType eps()          { return ScalarTraits<magnitudeType>::eps(); }
  static inline magnitudeType sfmin()        { return ScalarTraits<magnitudeType>::sfmin(); 
}
  static inline magnitudeType base()         { return ScalarTraits<magnitudeType>::base(); }
  static inline magnitudeType prec()         { return ScalarTraits<magnitudeType>::prec(); }
  static inline magnitudeType t()            { return ScalarTraits<magnitudeType>::t(); }
  static inline magnitudeType rnd()          { return ScalarTraits<magnitudeType>::rnd(); }
  static inline magnitudeType emin()         { return ScalarTraits<magnitudeType>::emin(); }
  static inline magnitudeType rmin()         { return ScalarTraits<magnitudeType>::rmin(); }
  static inline magnitudeType emax()         { return ScalarTraits<magnitudeType>::emax(); }
  static inline magnitudeType rmax()         { return ScalarTraits<magnitudeType>::rmax(); }
  static magnitudeType magnitude(ComplexT a)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        a, "Error, the input value to magnitude(...) a = " << a << " can not be NaN!" );
#endif
      return std::abs(a);
    }
  static inline ComplexT zero()              { return 
ComplexT(ScalarTraits<magnitudeType>::zero(),ScalarTraits<magnitudeType>::zero()); }
  static inline ComplexT one()               { return 
ComplexT(ScalarTraits<magnitudeType>::one(),ScalarTraits<magnitudeType>::zero()); }
  static inline ComplexT conjugate(ComplexT a){ return ComplexT(a.real(),-a.imag()); }
  static inline magnitudeType real(ComplexT a) { return a.real(); }
  static inline magnitudeType imag(ComplexT a) { return a.imag(); }
  static inline ComplexT nan()               { return 
ComplexT(ScalarTraits<magnitudeType>::nan(),ScalarTraits<magnitudeType>::nan()); }

  static inline bool isnaninf(ComplexT x)    { return 
ScalarTraits<magnitudeType>::isnaninf(x.real()) || 
ScalarTraits<magnitudeType>::isnaninf(x.imag()); }
  static inline void seedrandom(unsigned int s) { 
ScalarTraits<magnitudeType>::seedrandom(s); }
  static inline ComplexT random()
    {
      const <_MT_> rnd1 = ScalarTraits<magnitudeType>::random();
      const <_MT_> rnd2 = ScalarTraits<magnitudeType>::random();
      return ComplexT(rnd1,rnd2);
    }
  static inline std::string name() { return 
std::string("std::complex<")+std::string(ScalarTraits<magnitudeType>::name())+std::string(">"); 
}
  // This will only return one of the square roots of x, the other can be obtained by taking its conjugate
  static inline ComplexT squareroot(ComplexT x)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
        x, "Error, the input value to squareroot(...) x = " << x << " can not be NaN!" );
#endif
      typedef ScalarTraits<magnitudeType>  STMT;
      const <_MT_> r  = x.real(), i = x.imag(), zero = STMT::zero(), two = 2.0;
      const <_MT_> a  = STMT::squareroot((r*r)+(i*i));
      const <_MT_> nr = STMT::squareroot((a+r)/two);
      const <_MT_> ni = ( i == zero ? zero : STMT::squareroot((a-r)/two) );
      return ComplexT(nr,ni);
    }
  } // namespace Teuchos
#endif

#endif
