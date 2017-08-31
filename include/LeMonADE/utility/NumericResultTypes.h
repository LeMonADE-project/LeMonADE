/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_UTILITY_NUMERICRESULTTYPES_H
#define LEMONADE_UTILITY_NUMERICRESULTTYPES_H

/*****************************************************************************/
/**
 * @file
 * @brief class template NumericResultTypes, NumericTypeFlags, IfThenElse
 **/

/**
 * @file
 *
 * @namespace LemonadeHelper
 *
 * @brief The namespace <b>LeMonADeHelper</b> provides safety type-conversion-operations
 *
 * @details In this namespace all primitive type such Vector3D and safety-conversion-operator are defined
 *
 * @todo rename the namespace?
 * @todo rename the file? split it?
 **/
namespace LemonadeHelper
{


//! Generic template class providing the information numeric type is floating point with sign.
template < class T > struct NumericTypeFlags {};/*
{
    enum{
    is_float   = false,
    is_signed  = true
  };
};*/


//! Template for numeric type 'float'
template < >  struct NumericTypeFlags< float > { enum{ is_float   = true, is_signed  = true };};

//! Template for numeric type 'double'
template < >  struct NumericTypeFlags< double > { enum{  is_float   = true, is_signed  = true };};

//! Template for numeric type 'int8_t'
template < >  struct NumericTypeFlags < int8_t >   { enum{ is_float = false, is_signed   = true }; };

//! Template for numeric type 'int16_t'
template < >  struct NumericTypeFlags < int16_t >  { enum{ is_float = false, is_signed   = true }; };

//! Template for numeric type 'int32_t'
template < >  struct NumericTypeFlags < int32_t >  { enum{ is_float = false, is_signed   = true }; };

//! Template for numeric type 'int64_t'
template < >  struct NumericTypeFlags < int64_t >  { enum{ is_float = false, is_signed   = true }; };


//! Template for numeric type 'uint8_t'
template < >  struct NumericTypeFlags < uint8_t >   { enum{ is_float = false, is_signed   = false }; };

//! Template for numeric type 'uint16_t'
template < >  struct NumericTypeFlags < uint16_t >  { enum{ is_float = false, is_signed   = false }; };

//! Template for numeric type 'uint32_t'
template < >  struct NumericTypeFlags < uint32_t >  { enum{ is_float = false, is_signed   = false }; };

//! Template for numeric type 'uint64_t'
template < >  struct NumericTypeFlags < uint64_t >  { enum{ is_float = false, is_signed   = false }; };

//! Generic template class for if-then-clause for numeric types.
template < bool If, class Then, class Else > struct IfThenElse {};

//! Template for if-then-clause for 'then' branch.
template < class Then, class Else > struct IfThenElse < true , Then, Else > { typedef Then result;};

//! Template for if-then-clause for 'else' branch.
template < class Then, class Else > struct IfThenElse < false , Then, Else > { typedef Else result;};

//! Generic template class for providing which type has larger size.
template < class T1, class T2 > struct LargerType
{
   typedef typename IfThenElse < (sizeof (T1) > sizeof(T2)) , T1, T2 > ::result result;
};

//! Generic template class for providing the maximum between two uint.
template < uint v1, uint v2 > struct MaxValue
{
  enum { result = v1 > v2 ? v1 : v2 };
};

//! Generic template class for defining the result type of numeric operation.
template < class T1, class T2 = T1 > struct ResultTypes;

//! Generic template class for defining the return type of numeric operation.
template < class T > struct ResultTypes < T , T >
{
  typedef T stronger_type;
  typedef T product_type;
  typedef double double_type;

};

//! Generic template class providing the information of numeric type it's sign and number of bytes.
template < bool is_signed, uint size_of > struct SpecifyInt {};

//! Template for result type 'int8_t'
template <  > struct SpecifyInt < true, 1 > { typedef int8_t   result; };

//! Template for result type 'uint8_t'
template <  > struct SpecifyInt < false, 1 > { typedef uint8_t result; };

//! Template for result type 'int16_t'
template <  > struct SpecifyInt < true, 2 > { typedef int16_t   result; };

//! Template for result type 'uint16_t'
template <  > struct SpecifyInt < false, 2 > { typedef uint16_t result; };

//! Template for result type 'int32_t'
template <  > struct SpecifyInt < true, 4 > { typedef int32_t result; };

//! Template for result type 'uint32_t'
template <  > struct SpecifyInt < false, 4 > { typedef uint32_t result; };

//! Template for result type 'int64_t'
template <  > struct SpecifyInt < true, 8 > { typedef int64_t result; };

//! Template for result type 'uint64_t'
template <  > struct SpecifyInt < false, 8 > { typedef uint64_t result; };

}

namespace Lemonade
{
/**
 * @class NumericResultTypes
 * @brief Used by Vector3D to ensure well defined result types for operations involving different types
 **/
template < class T1, class T2 = T1  > struct NumericResultTypes
{
    typedef typename LemonadeHelper::IfThenElse
    <
	// if T1 is float
	LemonadeHelper::NumericTypeFlags <T1> ::is_float,

	// then
	  typename LemonadeHelper::IfThenElse
	  <
	    // if type 2 is float
	    LemonadeHelper::NumericTypeFlags <T2> ::is_float,
	      // then take the longer one
	      typename LemonadeHelper::LargerType<T1,T2>::result,
	      // else take the float T1
	      T1
	    // endif T2 is float
	    >::result,
	// else ( if T1 is not float )
	typename LemonadeHelper::IfThenElse
	<
	// if type T2 is float
	    LemonadeHelper::NumericTypeFlags <T2> ::is_float,

	    // then take T2
	    T2,

	    // else
	    // return a signed integer, if one of them is signed
	    // and at the same time the longest one
	    typename
	      LemonadeHelper::SpecifyInt
	      <
		LemonadeHelper::MaxValue < LemonadeHelper::NumericTypeFlags <T1> ::is_signed, LemonadeHelper::NumericTypeFlags <T2> ::is_signed > ::result,
		LemonadeHelper::MaxValue < sizeof(T1),sizeof(T2) > ::result
	      > ::result

	// endif T2 is float
	> :: result

    // endif T1 is float
    > :: result stronger_type;

    typedef stronger_type product_type;
    typedef double double_type;

};
}

#endif
