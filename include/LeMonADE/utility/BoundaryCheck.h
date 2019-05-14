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

#ifndef LEMONADE_UTILITY_BOUNDARYCHECK_H
#define LEMONADE_UTILITY_BOUNDARYCHECK_H

#include <stdexcept>
#include <iostream>
#include <sstream>

/**
 * @brief Index out of bounds exception to be used for arrays, extends the message by index and first element.
 *
 * @todo do we use this for DEBUG?
 **/
struct IndexOutOfBoundsException : public std::out_of_range
{
	std::string msg;
	template < class Array >
	IndexOutOfBoundsException(const Array& array, size_t idx) : std::out_of_range("Throwing std::out_of_range exception.")
	{
	  std::ostringstream strm;
	  strm << "Accessing array " << std::endl;
	  strm << "first element = " << array[0] << std::endl;
	  strm << "Index out of bounds, accessing array with key " << idx << " \n";
	  msg = strm.str();
	}

 	virtual ~IndexOutOfBoundsException() throw() {}
 	const char* what() const throw()
 	{
	  return std::string(msg + std::string(out_of_range::what())).c_str();
	}
};


/**
 * @brief A policy class to determine that boundaries are checked.
 *
 * @todo do we use this for DEBUG?
 **/
class CheckDynamicBounds
{
public:
  ///@brief Checks, if the array boundaries of "array" are respected by "key" and throws IndexOutOfBoundsException, if not.
  template < class Array, class KeyType >
  static void isValidIdx( const Array& array, const KeyType& key)
  {
	if ( size_t ( key ) >= array.size() )
	{
	  throw IndexOutOfBoundsException(array,key);
	}
  }
};

/**
 * @brief A policy class to determine that compile-time known boundaries are checked.
 *
 * @tparam size the compile-time known size of the array.
 * @tparam verbose If true, the "<<" output of the whole array is added to the exception message. Default: false.
 *
 * @todo do we use this for DEBUG?
 **/
template < uint size, bool verbose = false >
class CheckStaticBounds
{
public:
  //! Checks, if the array boundaries of "array" are respected by "key" and throws std::out_of_range, if not.
  template < class Array, class KeyType >
  static void isValidIdx( const Array& array, const KeyType& key)
  {
	if ( size_t ( key ) >= size )
	{
	  std::ostringstream msg;
	  if ( verbose ) msg << "container content: " << array << std::endl;
	  msg << "access with key " << key << std::endl;
	  throw std::out_of_range(msg.str().c_str());
	}
  }
};


/**
 * @brief A policy class to determine that not boundaries will be checked.
 *
 * @todo do we use this for DEBUG?
 */
class DontCheckBounds
{
public:
  ///@brief Empty function to fulfill the policy of having an isValidIdx() method. No boundaries are checked.
  template < class Array, class RefType >
  static void isValidIdx( const Array& a, const RefType& r){ }
};



#endif

