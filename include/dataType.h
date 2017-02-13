// Depending on the type 'T' in the getMpiDataType(T), this function will return back the correct MPI-data-type 
// This piece of code is mostly same as the following solution, with very few modifications.
// https://chuckaknight.wordpress.com/2013/03/13/intrinsic-type-conversion-using-template-specialization/

#ifndef __DATATYPE_H__
#define __DATATYPE_H__

#include <mpi.h>
#include <stdexcept>


#define getMpiDataType(T) ConvertType(getAbstractionDataType<T>())

namespace Abstraction
{
typedef enum 
	{
		type_unknown = 0,
		type_char,              //MPI_CHAR signed char
		type_unsigned_char,     //MPI_UNSIGNED_CHAR unsigned char
		type_short,             //MPI_SHORT signed short int
		type_unsigned_short,    //MPI_UNSIGNED_SHORT unsigned short int
		type_int,               //MPI_INT signed int
		type_unsigned_int,      //MPI_UNSIGNED unsigned int
		type_long,              //MPI_LONG signed long int
		type_unsigned_long,     //MPI_UNSIGNED_LONG unsigned long int
		type_float,             //MPI_FLOAT float
		type_double,             //MPI_DOUBLE double
    type_long_double        //MPI_LONG_DOUBLE long double
//		type_c_complex,     
//    type_c_double_complex
	}DataType;
}
 

template <class T>
Abstraction::DataType getAbstractionDataType()
{ throw std::runtime_error("Intrinsic type not supported by the abstraction."); }



template <>
inline Abstraction::DataType getAbstractionDataType<char>()
{ return Abstraction::type_char; }
								 
template <>
inline Abstraction::DataType getAbstractionDataType<unsigned char>()
{ return Abstraction::type_unsigned_char; }

template <>
inline Abstraction::DataType getAbstractionDataType<short>()
{ return Abstraction::type_short; }

template <>
inline Abstraction::DataType getAbstractionDataType<unsigned short>()
{ return Abstraction::type_unsigned_short; }

template <>
inline Abstraction::DataType getAbstractionDataType<int>()
{ return Abstraction::type_int;}
template <>
inline Abstraction::DataType getAbstractionDataType<unsigned int>()
{ return Abstraction::type_unsigned_int; }

template <>
inline Abstraction::DataType getAbstractionDataType<long>()
{ return Abstraction::type_long; }

template <>
inline Abstraction::DataType getAbstractionDataType<unsigned long>()
{ return Abstraction::type_unsigned_long; }

template <>
inline Abstraction::DataType getAbstractionDataType<float>()
{ return Abstraction::type_float; }

template <>
inline Abstraction::DataType getAbstractionDataType<double>()
{ return Abstraction::type_double; }

template <>
inline Abstraction::DataType getAbstractionDataType<long double>()
{ return Abstraction::type_long_double; }


//template <>
//inline Abstraction::DataType getAbstractionDataType<std::complex<int> >()
//{ return Abstraction::type_c_complex; }

//template <>
//inline Abstraction::DataType getAbstractionDataType<std::complex<double> >()
//{ return Abstraction::type_c_double_complex; }


template<class T>
void PrintTypeEnumerationValue()
{
	printf("Type enumeration value is %d\n", getAbstractionDataType<T>());
}

	
static MPI_Datatype ConvertType(Abstraction::DataType type)
{
  switch(type)
  {
    case Abstraction::type_char: return MPI_CHAR;
    case Abstraction::type_unsigned_char: return MPI_UNSIGNED_CHAR;
    case Abstraction::type_short: return MPI_SHORT;
    case Abstraction::type_unsigned_short: return MPI_UNSIGNED_SHORT;
    case Abstraction::type_int: return MPI_INT;
    case Abstraction::type_unsigned_int: return MPI_UNSIGNED;
    case Abstraction::type_long: return MPI_LONG;
    case Abstraction::type_unsigned_long: return MPI_UNSIGNED_LONG;
    case Abstraction::type_float: return MPI_FLOAT;
    case Abstraction::type_double: return MPI_DOUBLE;
    case Abstraction::type_long_double: return MPI_LONG_DOUBLE;
//    case Abstraction::type_c_complex: return MPI_C_COMPLEX;     // TODO: these complex types are not compatible with std::complex data types and require MPI-3.1+ implementation. Thus, checkpoint for std::complex data types has to be implemented separately 
//    case Abstraction::type_c_double_complex:{ std::cout << "returning MPI_C_DOUBLE_COMPLEX\n"; return MPI_C_DOUBLE_COMPLEX; }
  };
  	throw std::runtime_error("MPI_Datatype Convert(Abstraction::DataType) failed");
}

#endif
