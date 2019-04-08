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
 * @fn  MinImageDistanceComponentForPowerOfTwo
 * @brief calculates the minimal distances of images for one component 
 * @return int 
 * @param x1 absolute coordinate
 * @param x2 absolute coordinate
 * @param LatticeSize size of the box in the direction of the given coordinates
 */
inline int MinImageDistanceComponentForPowerOfTwo(const int x1, const int x2, const uint32_t latticeSize )
{
	//this is only valid for absolute coordinates
	uint32_t latticeSizeM1(latticeSize-1);
	return ( (((x2-x1)&latticeSizeM1) < (latticeSize/2)) ? ((x2-x1) & latticeSizeM1) :  -((x1-x2) & latticeSizeM1));
} 
/**
 * @fn MinImageVectorForPowerOfTwo
 * @brief returns the shortest vector between two images for a box size of power of 2 
 * @return vector 
 * @param R1 position vector  
 * @param R2 position vector 
 * @param ing container containing (all) system information
 */
template < class IngredientsType>
VectorInt3 MinImageVectorForPowerOfTwo (const VectorInt3 R1, const VectorInt3 R2, IngredientsType& ing)
{
  VectorInt3 dist;
  dist.setX(MinImageDistanceComponentForPowerOfTwo(R1.getX(),R2.getX(),ing.getBoxX()));
  dist.setY(MinImageDistanceComponentForPowerOfTwo(R1.getY(),R2.getY(),ing.getBoxY()));
  dist.setZ(MinImageDistanceComponentForPowerOfTwo(R1.getZ(),R2.getZ(),ing.getBoxY()));
  return dist;
}
/**
 * @fn MinImageDistanceForPowerOfTwo
 * @brief returns the minimal distance between two images for a box size of power of 2 
 * @return double 
 * @param R1 position vector  
 * @param R2 position vector 
 * @param ing container containing (all) system information
 */
template < class IngredientsType>
double MinImageDistanceForPowerOfTwo (const VectorInt3 R1, const VectorInt3 R2, IngredientsType& ing)
{
  return MinImageVectorForPowerOfTwo(R1,R2,ing).getLength();
}

// Implementation for arbitrary box dimensions. 

/**
 * @fn MinImageDistance
 * @brief returns the minimal distance between two images for an arbitrary box size 
 * @return double 
 * @param R1 position vector  
 * @param R2 position vector 
 * @param ing container containing (all) system information
 */
template < class IngredientsType>
double MinImageDistance (const VectorInt3 R1, const VectorInt3 R2, IngredientsType& ing)
{
  return MinImageVector(R1,R2,ing).getLength();
}
  
/**
 * @fn MinImageVector
 * @brief returns the shortest vector between two images for an arbitrary box size 
 * @return vector 
 * @param R1 position vector  
 * @param R2 position vector 
 * @param ing container containing (all) system information
 */
template < class IngredientsType>
VectorInt3 MinImageVector (const VectorInt3 R1, const VectorInt3 R2, IngredientsType& ing)
{
  VectorInt3 dist;
  dist.setX(MinImageDistanceComponen(R1.getX(),R2.getX(),ing.getBoxX()));
  dist.setY(MinImageDistanceComponen(R1.getY(),R2.getY(),ing.getBoxY()));
  dist.setZ(MinImageDistanceComponen(R1.getZ(),R2.getZ(),ing.getBoxY()));
  return dist;
}
/**
 * @fn  MinImageDistanceComponen
 * @brief calculates the minimal distances of images for an arbitrary box size 
 * @return int 
 * @param x1 absolute coordinate
 * @param x2 absolute coordinate
 * @param LatticeSize size of the box in the direction of the given coordinates
 */
inline int MinImageDistanceComponen(const int x1, const int x2, const uint32_t latticeSize )
{
	//this is only valid for absolute coordinates
	int distance(x2-x1);
	int latticeHalf(latticeSize/2);
// 	while( std::abs(distance) > latticeHalf ) { distance=distance +- std::signbit() latticeSize; }
	if      (  distance > latticeHalf ) { while(-(distance-latticeSize) < latticeHalf ) distance -= latticeSize; }
	else if ( -distance > latticeHalf ) { while( (distance+latticeSize) < latticeHalf ) distance += latticeSize; }
	return distance;
}
};


#endif
