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

#ifndef LEMONADE_UTILITY_FASTBONDSET_H
#define LEMONADE_UTILITY_FASTBONDSET_H

#include <map>
#include <iostream>
#include <utility>
#include <sstream>

#include <LeMonADE/utility/Vector3D.h>

/***********************************************************/
/**
 * @file
 * @class FastBondset
 *
 * @brief Contains and manages the used set of bond-vectors.
 *
 * @details This FastBondset class allows bond vectors with components \a x , \a y ,\a z
 * to be -4 <= x,y,z <=4 for fast look-up. For larger allowed bond vectors use
 * class SlowBondset and change FeatureBondeset, appropriately.
 * Also, to improve performance, some hardly used safety checks were omitted
 * here, which are present in SlowBondset.
 *
 */
/***********************************************************/

class FastBondset
{
private:

	//! Look-up table telling if a certain bond-vector is valid.
	bool bondsetLookup[512];

	//! Translates bond-vector to index in the (fast) look-up table (bondsetLookup).
	uint32_t bondVectorToIndex(const VectorInt3& bondVector) const;

protected:

	//! Each bond-vector is represented by a std::map <int32_t, VectorInt3>, where the identifier is stored as int32_t
	std::map < int32_t , VectorInt3> BondVectors;

	//! True if bondsetLookup is already synchronized with added allowed bond-vectors, false otherwise.
	bool lookupSynchronized;

public:

	//! Default constructor
	FastBondset();

	//! Copy constructor
	FastBondset(const FastBondset& copyBondSet);

	//! Default destructor
	virtual ~FastBondset();

	//! Standard iterator for map of BondVectors
	typedef std::map <int32_t, VectorInt3>::const_iterator iterator;

	//! Get identifier by the given bond.
	int32_t getBondIdentifier(int32_t x,int32_t y,int32_t z) const;

	//! Get bond by its identifier (character)
	const VectorInt3 getBondVector(int32_t identifier) const;

	//! Adding bond to the allowed set
	void addBond(VectorInt3 bondVector, int32_t identifier);

	//! Adding bond to the allowed set
	void addBond(int32_t x, int32_t y, int32_t z, int32_t identifier);


	//! Updates the look-up table of bond-vectors
	void updateLookupTable();

	//! Resets/delete the look-up table of bond-vectors
	void resetLookupTable();

	//! Check if a vector is a valid bond-vector (i.e. part of the set)
	bool isValid(const VectorInt3& bondVector) const;

	//! Check if a vector is a valid bond-vector (i.e. part of the set)
	bool isValidStrongCheck(const VectorInt3& bondVector) const;


	//! Clear the look-up table and storing map of bond-vectors
	void clear();


	//! Fct for providing iteration begin of the map from outside
	iterator begin() const {return BondVectors.begin();}

	//! Fct for providing iteration end of the map from outside
	iterator end() const {return BondVectors.end();}

	//! Getting the size of the map from outside
	size_t size() const {return BondVectors.size();}

	//! Adding the predefined set of \b bccBFM-bond-vectors to the map
	void addBccBFMclassicBondset();

	//! Adding the predefined set of \b scBFM-bond-vectors to the map
	void addBFMclassicBondset();

};

//////////////////////////////////////////////////////////////////////////////
// inline functions defined here
//////////////////////////////////////////////////////////////////////////////


/**
 * @details There's is no boundary-check esp -3 <= x,y,z <= 3. This is maybe critical
 *          if you want to check if a bond is allowed
 *
 * @param bondVector Reference to VectorInt3 as bond-vector to check.
 * @return True if bond-vector is allowed, false otherwise.
 */
inline bool FastBondset::isValid(const VectorInt3& bondVector) const
{
	return bondsetLookup[bondVectorToIndex(bondVector)];
}

/**
 * @details There's is boundary-check makes sure that -3 <= x,y,z <= 3. This is useful
 *          if you want to check if a bond is allowed
 *
 * @param bondVector Reference to VectorInt3 as bond-vector to check.
 * @return True if bond-vector is allowed, false otherwise.
 */
inline bool FastBondset::isValidStrongCheck(const VectorInt3& bondVector) const
{
	if(bondVector.getX() < -3)
		return false;

	if(bondVector.getX() > 3)
		return false;

	if(bondVector.getY() < -3)
		return false;

	if(bondVector.getY() > 3)
		return false;

	if(bondVector.getZ() < -3)
		return false;

	if(bondVector.getZ() > 3)
		return false;




	return bondsetLookup[bondVectorToIndex(bondVector)];
}


/**
 * @details Translates a bond-vector into the corresponding lookup table index.
 * The lookup method is based on the fact that for integers in the range
 * -4<= int <=4 there is a well defined map between the first three bits of the integers
 * and the number( except for 4,(-4), but this is always rejected):
 * * -4 &rarr; (-4 &7) = 4
 * * -3 &rarr; (-3 &7) = 5
 * * -2 &rarr; (-2 &7) = 6
 * * -1 &rarr; (-1 &7) = 7
 * *  0 &rarr; ( 0 &7) = 0
 * *  1 &rarr; ( 1 &7) = 1
 * *  2 &rarr; ( 2 &7) = 2
 * *  3 &rarr; ( 3 &7) = 3
 * *  4 &rarr; ( 4 &7) = 4
 *
 * @param bondVector Reference to VectorInt3 as bond-vector to put on fast look-up-table.
 * @return A unique number representing the bond-vector.
 */
inline uint32_t FastBondset::bondVectorToIndex(const VectorInt3& bondVector) const
{
	return (bondVector.getX() & 7) + ((bondVector.getY() &7) << 3) + ((bondVector.getZ() &7) << 6);
}

#endif /* LEMONADE_UTILITY_FASTBONDSET_H */
