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

#ifndef LEMONADE_UTILITY_DISTANCECALCULATION_H
#define LEMONADE_UTILITY_DISTANCECALCULATION_H

#include <iostream>

#include "extern/loki/NullType.h"

#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/SafeCast.h>

namespace Lemonade
{

/**
 * @brief Calculation of Distance in the Minimum Image Convention.
 *
 * @deprecated
 *
 * @param[out] distance Distance in MIC
 * @param[in] period Length of the Unit Cell
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Different namespace?
 **/
template < class T1, class T2 > void reduceDistanceInPeriodicSpace( T1& distance, const T2& period)
{
	if (  distance > period/2 ) {while(-(distance-period)<period/2) distance -= Lemonade::safe_cast<T1>(period);  return; }
	else
	if ( -distance > period/2 ) { while((distance+period)<period/2) distance += Lemonade::safe_cast<T1>(period); return; }
}


/**
 * @brief Calculation of Vector in Cartesian space using the Minimum Image Convention.
 *
 * @param a Vector in Cartesian space
 * @param b Vector in Cartesian space
 * @param box Length of the Unit Cell
 * @return Vector under MIC
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Different namespace?
 **/
template < class T, class BoxType >
Vector3D<T> calcDistanceVector3D (const Vector3D< T >& a, const Vector3D< T >& b, const BoxType& box )
{
	Vector3D < T > result (b-a);
	if ( box.isPeriodicX() ){ reduceDistanceInPeriodicSpace ( result[0], box.getBoxX() ); }
	if ( box.isPeriodicY() ){ reduceDistanceInPeriodicSpace ( result[1], box.getBoxY() ); }
	if ( box.isPeriodicZ() ){ reduceDistanceInPeriodicSpace ( result[2], box.getBoxZ() ); }
	return result;
}


/**
 * @brief Calculation of difference vector in Cartesian space.
 *
 * @param a Vector in Cartesian space
 * @param b Vector in Cartesian space
 * @param box Length of the Unit Cell
 * @return Difference vector
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Different namespace?
 **/
template < class T >
Vector3D<T> calcDistanceVector3D (const Vector3D< T >& a, const Vector3D< T >& b, const Loki::NullType& box )
{
	Vector3D < T > result (b-a);
	return result;
}

};


#endif
