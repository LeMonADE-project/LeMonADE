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

#ifndef LEMONADE_UTILITY_SAFECAST_H
#define LEMONADE_UTILITY_SAFECAST_H

/*****************************************************************************/
/**
 * @file
 * @brief Definitions of class templates SafeCastBackend and safe_cast
 * */
/*****************************************************************************/

#include <stdexcept>

/*****************************************************************************/
/**
 * @class SafeCastBackend
 *
 * @brief Backend for safe_cast, making sure some casts are forbidden
 **/
namespace LemonadeHelper
{

//! Generic template class providing a static_cast to \a From to \a To
template < class From, class To > struct SafeCastBackend
{
  static To Cast( From& src) { return static_cast<To>(src); };
};

//! Generic template class providing a static_cast to \a To to \a To
template < class To > struct SafeCastBackend<To,To>
{
  static const To& Cast( const To& src) { return src;};
};

/// specializations explicitly forbidding some casts:

// Template for explicitly forbidding some casts

//! Template for explicitly forbidding cast from \a int64_t to \a int
template <> struct SafeCastBackend<int64_t,int> {};

//! Template for explicitly forbidding cast from \a int64_t to \a uint
template <> struct SafeCastBackend<uint64_t,uint> {};

//! Template for explicitly forbidding cast from \a double to \a float
template <> struct SafeCastBackend<double,float> {};

//! Template for explicitly forbidding cast from \a float to \a int
template <> struct SafeCastBackend<float,int> {};

//! Template for explicitly forbidding cast from \a double to \a int
template <> struct SafeCastBackend<double,int> {};

//! Template for explicitly forbidding cast from \a float to \a uint
template <> struct SafeCastBackend<float,uint> {};

//! Template for explicitly forbidding cast from \a double to \a uint
template <> struct SafeCastBackend<double,uint> {};
};

namespace Lemonade
{
  /*****************************************************************************/
  /**
  * @class safe_cast
  *
  * @brief Provides save type-cast, i.e. only without loss of information
  **/
  template < class To > struct safe_cast
  {
    To tmp;

    template < class From > safe_cast(From  src):tmp( LemonadeHelper::SafeCastBackend< From, To >::Cast( src ) ){}

    operator To () const {return tmp;}
  };
};

#endif
