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

#ifndef LEMONADE_UTILITY_VECTOR3D_H
#define LEMONADE_UTILITY_VECTOR3D_H

/**
 * @file Vector3D.h
 * @brief  Template definition for a general 3D vector
 *
 * @author werner
 * @date   12.03.2013
 **/
/****************************************************************/

#include <iostream>
#include <cmath>
#include <stdint.h>
#include <string>
#include <sys/types.h>
#include <iomanip>

#include <LeMonADE/utility/BoundaryCheck.h>
//#include <LeMonADE/utility/StaticArray.h>
#include <LeMonADE/utility/NumericResultTypes.h>
#include <LeMonADE/utility/SafeCast.h>

namespace LemonadeHelper
{

template < class T >
struct SafeCastPolicy
{
  template < class T1 > 
  static T cast(T1 val) { return Lemonade::safe_cast< T >(val); }
};

}

using Lemonade::NumericResultTypes;
using Lemonade::safe_cast;

// namespace Lemonade
// {

//@details The Vector3D contains three coordinates x, y and z and one
// * component, w, as metadata. Vector3D<T> bundles simple arithmetic operations
// * affecting three coordinates x,y and z. The fourth component w can be used for
// * meta data and leads to an improved layout in memory especially for vectors
// * of Vector3D<T>s .

/**
 * @class Vector3D
 *
 * @brief \b An \b Euclidean-Vector \b implementation \b in \b Cartesian-space \b (3D) \b of \b arbitrary \b primitive \b numeric \b type.
 *
 * @details The Vector3D contains three coordinates x, y and z.
 * Vector3D<T> bundles simple arithmetic operations affecting three coordinates x,y and z.
 *
 *
 * @todo theres not metadata w anymore -but old documentation?
 * @tparam T built-in-type forming the x, y an z component.
**/
template <class T>
class Vector3D
{
private:
  
	//! Providing a safe cast between Vector3D of different types
	typedef LemonadeHelper::SafeCastPolicy<T> Casting;
	
#ifdef DEBUG
	typedef CheckStaticBounds<3> BoundaryCheck;
#else
	typedef DontCheckBounds BoundaryCheck;
#endif

	/**
	 * @union
	 * @brief Storage of the three Cartesian-components
	 */
	union
	{
		struct
		{
			T X; //!< x-coordinate in the Cartesian space
			T Y; //!< y-coordinate in the Cartesian space
			T Z; //!< z-coordinate in the Cartesian space
		};
 		T v[3]; //!< x,y,z-coordinate in the Cartesian space as array
	};
	
	/**
	 * @brief Read-only access of the \a i-th element in Vector3D
	 *
	 * @param i i=0:x-coordinate \n i=1:y-coordinate \n i=2:z-coordinate
	 * @return Value at specified coordinate
	 */
	const T& const_access (uint i) const
	{
	  BoundaryCheck::isValidIdx(*this,i); 
	  return v[i];
	}
	
	/**
     * @brief Read-Write access of the \a i-th element in Vector3D
	 *
	 * @param i i=0:x-coordinate \n i=1:y-coordinate \n i=2:z-coordinate
	 * @return Value at specified coordinate
	 */
	T& access (uint i) 
	{
	  BoundaryCheck::isValidIdx(*this,i); 
	  return v[i];
	}
	

public:
  
// 	typedef T coordinate_type;

	/// get the ith component i=0 --> X, ... .

	/**
	 * @brief Read-only access of the \a i-th element in Vector3D
	 *
	 * @param i i=0:x-coordinate \n i=1:y-coordinate \n i=2:z-coordinate
	 * @return Value at specified coordinate
	 **/
 	const T& getCoordinate(unsigned i)	const {return this->const_access(i);}

	/**
	 * @brief Read-only access of the \a x-coordinate in Vector3D
	 *
	 * @return Value at x-coordinate
	 **/
	const T& getX() 			const {return X;}

	/**
	 * @brief Read-only access of the \a y-coordinate in Vector3D
	 *
	 * @return Value at y-coordinate
	 **/
	const T& getY() 			const {return Y;}

	/**
	 * @brief Read-only access of the \a z-coordinate in Vector3D
	 *
	 * @return Value at z-coordinate
	 **/
	const T& getZ() 			const {return Z;}



	/// Analog to get(i).
	/**
     * @brief Write access to set \a val of the \a i-th element in Vector3D
	 *
	 * @param i i=0:x-coordinate \n i=1:y-coordinate \n i=2:z-coordinate
	 * @param val Value to set at specified coordinate
	 **/
	void setCoordinate(unsigned i, const T& val){ this->access(i) = val; }
	
	/**
     * @brief Write access to set all element in Vector3D
	 *
	 * @param i i=0:x-coordinate \n i=1:y-coordinate \n i=2:z-coordinate
	 * @param val Value to set at specified coordinate
	 **/
	template < class T1, class T2, class T3 >
	void setAllCoordinates(const T1& x, const T2& y, const T3& z)
	{
	  X = Casting::cast(x); 
	  Y = Casting::cast(y); 
	  Z = Casting::cast(z);
	}
	
	/**
	 * @brief Write access to set \a val on the \a x-coordinate in Vector3D
	 *
	 * @param val Value to set at x-coordinate
	 **/
	template < class T1 > void setX(const T1& val) {X = Casting::cast(val);}

	/**
	 * @brief Write access to set \a val on the \a y-coordinate in Vector3D
	 *
	 * @param val Value to set at y-coordinate
	 **/
	template < class T1 > void setY(const T1& val) {Y = Casting::cast(val);}

	/**
	 * @brief Write access to set \a val on the \a z-coordinate in Vector3D
	 *
	 * @param val Value to set at z-coordinate
	 **/
	template < class T1 > void setZ(const T1& val) {Z = Casting::cast(val);}
	
	//! Components can be accessed directly by using Vector3D < T > as an array.
	const T& operator [](unsigned i)const{ return this->const_access(i);}

	//! The non-const version of operator[], probably it will/should be removed again.
	T& operator [](unsigned i){return  this->access(i);}

	//! Standard constructor. Default: all elements are set to Zero.
	Vector3D():X(0),Y(0),Z(0){}

	/**
	 * @brief Single-parameter copy-constructor for arbitrary Vector3D-like classes.
	 *
	 * @details The results of src.getX(), src.getY(), src.getZ() are type-casted and used to initialize the own components.
	 * @param src
	 */
	template < class VectorType > Vector3D(const VectorType & src)
	:X( Casting::cast(src.getX()) )
	,Y( Casting::cast(src.getY()) )
	,Z( Casting::cast(src.getZ()) )
	{}
	
	/**
	 * @brief Assignment operator providing deep-copy of the elements of /a src  to Vector3D.
	 *
	 * @details  The results of src.getX(), src.getY(), src.getZ() are type-casted and used to initialize the own components.
	 * @param src Vector3D to assign from
	 * @return Vector3D<T1> with deep-copy of \a src.
	 */
	template < class T1 > Vector3D& operator = ( const Vector3D < T1 > &src )
	{
		setX(Casting::cast(src.getX()));
		setY(Casting::cast(src.getY()));
		setZ(Casting::cast(src.getZ()));
		return *this;
	}
	
	// 3D constructor.
	/**
	 * @brief Standard constructor. Initializing all coordinates/elements of Vector3D with given parameter \a x,y,z.
	 *
	 * @details  The value of \a x,y,z are type-casted and used to initialize the own components.
	 *
	 * @param x x-coordinate in the Cartesian space
	 * @param y y-coordinate in the Cartesian space
	 * @param z z-coordinate in the Cartesian space
	 **/
	template < class T1, class T2, class T3 > Vector3D ( T1 x, T2 y, T3 z)
	:X(Casting::cast(x))
	,Y(Casting::cast(y))
	,Z(Casting::cast(z))
	{}

	/**
	 * @brief \b Dot \b (Scalar) \b product returning a sum of type T.
	 *
	 * @details Components from \a src are safety type-casted to T before multiplying.
	 *
	 * @param src object of class Vector3D < src_T >
	 *
	 * @return Sum of the scalar product.
	 **/
	template < class T1 > 
	typename NumericResultTypes < T, T1 >::product_type operator * (const Vector3D < T1 > & src) const
	{
		typedef typename NumericResultTypes <T, T1 >::product_type ResultType;
		ResultType sum = 0;
		sum   = safe_cast<ResultType>(X) * safe_cast<ResultType>(src.getX());
		sum  += safe_cast<ResultType>(Y) * safe_cast<ResultType>(src.getY());
		sum  += safe_cast<ResultType>(Z) * safe_cast<ResultType>(src.getZ());
		
		return sum;
	}


	/**
	 * @brief \b Scalar \b multiplication of the Vector3D with \a multiplier.
	 *
	 * @details Returns a newly constructed Vector3D object as copy of the left hand side,
	 * which is multiplied by \a multiplier.
	 *
	 * @param multiplier Scalar to multiply with Vector3D
	 * @return Scaled Vector3D*multiplier
	 **/
	template < class T1 > Vector3D operator * (const T1& multiplier) const
	{
		Vector3D result(*this);
		result *= Casting::cast(multiplier);
		return result;
	}
	
	/**
	 * @brief \b Scalar \b multiplication ('division') of the Vector3D with \a 1.0/divisor.
	 *
	 * @details Returns a newly constructed Vector3D object as copy of the left hand side,
	 * which is multiplied by \a 1.0/divisor.
	 *
	 * @param divisor Scalar to divide the Vector3D
	 * @return Scaled Vector3D/divisor
	 *
	 * @todo can throw DividedByZero-Exception
	 **/
	template < class T1 > Vector3D operator / (const T1& divisor) const
	{
		return (*this)*(1.0/divisor);
	}
	
	/**
	 * @brief \b Subtraction of the Vector3D \a src from this Vector3D.
	 *
	 * @details Returns a newly constructed Vector3D object as copy of the left hand side,
	 * which is subtracted by elements of \a src.
	 *
	 * @param src object of class Vector3D < T1 >
	 * @return Subtracted Vector3D = this-\a src with stronger_type.
	 *
	 **/
	template < class T1 > 
	typename NumericResultTypes < Vector3D, Vector3D< T1 > >::stronger_type operator - (const Vector3D<T1>& src) const
	{
	    ///@todo prevent overflow for unsigned types
	    return typename NumericResultTypes < Vector3D, Vector3D< T1 > >::stronger_type (
	       getX()-src.getX()
	      ,getY()-src.getY()
	      ,getZ()-src.getZ()
	    );
	}
	
	/**
	 * @brief \b Addition of the Vector3D \a src to this Vector3D.
	 *
	 * @details Returns a newly constructed Vector3D object as copy of the left hand side,
	 * which is added by elements of \a src.
	 *
	 * @param src object of class Vector3D < T1 >
	 * @return Subtracted Vector3D = this+\a src with stronger_type.
	 *
	 **/
	template < class T1 > 
	typename NumericResultTypes < Vector3D, Vector3D< T1 > >::stronger_type operator + (const Vector3D<T1>& src) const
	{
	    return typename NumericResultTypes < Vector3D, Vector3D< T1 > >::stronger_type (
	       getX()+src.getX()
	      ,getY()+src.getY()
	      ,getZ()+src.getZ()
	    );
	}

	/**
	 * @brief \b Subtraction of the Vector3D \a src from this Vector3D.
	 *
	 * @details The value of elements of \a src are subtracted from \a x,y,z.
	 *
	 * @param src object of class Vector3D
	 *
	 **/
	void operator -= (const Vector3D  & src)
	{
		X -= src.getX();
		Y -= src.getY();
		Z -= src.getZ();
	}

	/*bool operator > (const Vector3D& rhs) const
	{
		T sqr_lhs((*this)*(*this)), sqr_rhs(rhs*rhs);
		return sqr_lhs > sqr_rhs;
	}

	bool operator < (const Vector3D& rhs) const
	{
		T sqr_lhs((*this)*(*this)), sqr_rhs(rhs*rhs);
		return sqr_lhs < sqr_rhs;
	}*/

	/**
	 * @brief \b Addition of the Vector3D \a src from this Vector3D.
	 *
	 * @details The value of elements of \a src are added to \a x,y,z.
	 *
	 * @param src object of class Vector3D
	 *
	 **/
	void operator += (const Vector3D& rhs)
	{
		X += rhs.getX();
		Y += rhs.getY();
		Z += rhs.getZ();
	}

	/**
	 * @brief \b Scalar \b multiplication of this Vector3D with \a multiplier.
	 *
	 * @details The value of these \a x,y,z are multiplied by a type-casted \a multiplier.
	 *
	 * @param multiplier Scalar to multiply with Vector3D
	 **/
	template < class T1 > Vector3D& operator *= (const T1& multiplier)
	{
		X = X * Casting::cast( multiplier );
		Y = Y * Casting::cast( multiplier );
		Z = Z * Casting::cast( multiplier );
		return *this;
	}

	/**
	 * @brief \b Scalar \b multiplication ('division') of this Vector3D with \a 1.0/divisor.
	 *
	 * @details The value of these \a x,y,z are multiplied by a type-casted \a 1.0/divisor.
	 *
	 * @param divisor Scalar to divide the Vector3D
	 **/
	template < class T1 > Vector3D& operator /= (const T1& divisor)
	{
		X = X / Casting::cast(divisor);
		Y = Y / Casting::cast(divisor);
		Z = Z / Casting::cast(divisor);
		return *this;
	}

	/**
	 * @brief \b Equality \b operator for Vector3D comparing the equality of all elements.
	 *
	 * @param src object of class Vector3D
	 * @return True if all elements are the equal. False, otherwise.
	 **/
	bool operator == (const Vector3D& src) const
	{
		if (getX() != src.getX()) return false;
		if (getY() != src.getY()) return false;
		if (getZ() != src.getZ()) return false;
		return true;
	}

	/**
	 * @brief \b Inequality \b operator for Vector3D comparing the equality of all elements.
	 *
	 * @param src object of class Vector3D
	 * @return True if at least on element is different. False, otherwise.
	 **/
	bool 	operator != (const Vector3D& src) const
	{
		return !(*this == src);
	}

	/**
	 * @brief Returns the Euclidean norm of the Vector3D as primitive numeric type double.
	 *
	 * @details Internally it's calculates its square with its primitive numeric type T
	 * and converts the results with cmath::sqrt to its norm.
	 *
	 * @return Euclidean norm L=(x*x+y*y+z*z)^0.5
	 *
	 * @todo There can be happen an overflow of the types in the multiplication.
	 */
	double getLength() const
	{
		T sqr_norm((*this)*(*this));
		return 	sqrt((double)sqr_norm);
	}

	/**
	 * @brief Normalization of this Vector3D to an Unit Vector of length equals Unity.
	 *
	 * @details Internally it's calculates its norm with its primitive numeric type T
	 * and converts the results with cmath::sqrt and divides itself by the result (length).
	 *
	 *
	 * @todo There can be happen an overflow of the types in the multiplication.
	 */
	Vector3D& normalize()
	{
	  T prod = (*this)*(*this); 
	  *this /= sqrt(prod);
	  if ( this->getX() != this->getX() || this->getY() != this->getY() || this->getZ() != this->getZ() )
	  {
		  std::ostringstream strm;
		  strm << "Vector3D::normalize() leading to nan " << std::endl;
		  throw std::runtime_error(strm.str().c_str());
	  }
	  return *this;
	}

	//! Returns the const Vector3D itself
	const Vector3D& getVector3D() const 
	{
		  return *this;
	}
	
	//! Returns the modifiable Vector3D itself
	Vector3D& modifyVector3D() 
	{
		  return *this;
	}

};

/**
 * @brief \b Cross \b product of the two Vector3D a %times; b.
 *
 * @details Returns a newly constructed Vector3D object as result of the cross product
 *
 * @param a object of class Vector3D
 * @param b object of class Vector3D
 * @return cross product a %times; b.
 *
 * @throw <std::runtime_error> if NAN happens
 *
 * @todo implement the cross-product with %-operator
 **/
template <class T>
T crossProduct(const T & a, const T b)
{
	T	 prod(	a.getY()*b.getZ() - a.getZ()*b.getY(),
			a.getZ()*b.getX() - a.getX()*b.getZ(),
			a.getX()*b.getY() - a.getY()*b.getX() );
	if ( prod.getX() != prod.getX() || prod.getY() != prod.getY() || prod.getZ() != prod.getZ() )
	{
	      std::ostringstream strm;
	      strm << "vectors leading to nan: (" << a << ") (" << b << ")" << std::endl;
	      throw std::runtime_error(strm.str().c_str());
	}
	else return prod;
}

/**
 * @brief \b Stream \b Out \b operator of the Vector3D
 *
 * @details Streams out the elements x,y,z of the vector separated by Space
 *
 * @param stream output-stream
 * @param p object of class Vector3D
 * @return output-stream
 **/
template < class T >
std::ostream& operator << (std::ostream& stream, const Vector3D < T > & p)
{
	stream 
//  	<< std::scientific
//  	<< std::setprecision(15)
	<<  p.getX() << " " <<  p.getY() << " " << p.getZ();
	return stream;
}

/**
 * @brief \b Stream \b In \b operator of the Vector3D
 *
 * @details Streams in the elements. And setting at first x, then y, at last z.
 *
 * @param stream input-stream
 * @param p object of class Vector3D
 * @return input-stream
 **/
template < class T >
std::istream& operator >> (std::istream& stream, Vector3D < T > & p)
{
	T temp;
	stream >> temp; p.setX(temp);
	stream >> temp; p.setY(temp);
	stream >> temp; p.setZ(temp);
	return stream;
}

/**
 * @brief \b Scalar \b multiplication of given Vector3D with \a scalar.
 *
 * @details The value of these \a x,y,z are multiplied by a type-casted \a scalar.
 *
 * @param scalar Scalar to multiply with Vector3D
 * @param vec object of class Vector3D
 * @todo type-casted scalar?
 **/
template < class T > Vector3D<T> operator * ( const T& scalar, const Vector3D<T>& vec )
{
    Vector3D<T> result(vec); result *= scalar; return result;
}

// template < class T1, class T2 > 
// Vector3D < typename NumericResultTypes < T1, T2 >::product_type > 
// operator * (const T1& scalar, const Vector3D<T2>& vec )
// {
//   typedef Vector3D< typename NumericResultTypes < T1, T2 >::product_type > result (vec);
//   result *= scalar;
//   return result;
// };

// Typedefs:

// most important:
//! Typedef of VectorInt3 as Vector3D < int32_t >
typedef Vector3D < int32_t > VectorInt3;

//! Typedef of VectorUint3 as Vector3D < uint32_t >
typedef Vector3D < uint32_t > VectorUint3;

//! Typedef of VectorFloat3 as Vector3D < float >
typedef Vector3D < float > VectorFloat3;

//! Typedef of VectorDouble3 as Vector3D < double >
typedef Vector3D < double > VectorDouble3;

//! Typedef of VectorLong3 as Vector3D < int64_t >
typedef Vector3D < int64_t > VectorLong3;

//! Typedef of VectorUlong3 as Vector3D < uint64_t >
typedef Vector3D < uint64_t > VectorUlong3;

///@todo clarify, if the following is confusing:

//! Typedef of VectorChar3 as Vector3D < int8_t >
typedef Vector3D < int8_t > VectorChar3;

//! Typedef of VectorUchar3 as Vector3D < uint8_t >
typedef Vector3D < uint8_t > VectorUchar3;

//! Typedef of VectorShort3 as Vector3D < int16_t >
typedef Vector3D < int16_t > VectorShort3;

//! Typedef of VectorUshort3 as Vector3D < uint16_t >
typedef Vector3D < uint16_t > VectorUshort3;



// template < class T > struct NumericTypeTraits < Vector3D < T > > { };
namespace Lemonade
{

/**
 * @todo merge with file NumericResultTypes?
 */
template < class T1, class T2 > struct NumericResultTypes < Vector3D <T1>, Vector3D <T2> >
{
  typedef Vector3D < typename NumericResultTypes<T1,T2>::stronger_type > stronger_type;
  typedef typename NumericResultTypes<T1,T2>::stronger_type product_type;
  typedef Vector3D < typename NumericResultTypes< stronger_type, stronger_type >::double_type > double_type;
};

template < class T > struct NumericResultTypes < Vector3D <T>, Vector3D <T> >
{
  typedef Vector3D < T > stronger_type;
  typedef T product_type;
  typedef Vector3D < typename NumericResultTypes< T >::double_type > double_type;
};

}
// }


#endif /* LEMONADE_UTILITY_VECTOR3D_H */

