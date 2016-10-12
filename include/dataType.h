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
		type_char,
		type_unsigned_char,
		type_short,
		type_unsigned_short,
		type_int,
		type_unsigned_int,
		type_long,
		type_unsigned_long,
		type_float,
		type_double
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
  };
  //	throw std::runtime_error("MPI_Datatype Convert(Abstraction::DataType) failed");
}

#endif
