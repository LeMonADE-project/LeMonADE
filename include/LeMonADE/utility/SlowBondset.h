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

#ifndef LEMONADE_UTILITY_SLOWBONDSET_H
#define LEMONADE_UTILITY_SLOWBONDSET_H

#include <map>
#include <iostream>
#include <utility>
#include <sstream>


#include <LeMonADE/utility/FastBondset.h>
#include <LeMonADE/utility/Vector3D.h>


/***********************************************************/
/**
 * @file
 * @class SlowBondset
 * @brief Contains and manages the used set of bondvectors.
 *
 * @details Contains and manages the used set of bondvectors. This class can
 * be used instead of the class Bondset, if bonds with components larger
 * than 4 or smaller than -4 are needed. It is more generally applicable
 * than Bondset, but is slower when checking bonds for validity.
 *
 */
/***********************************************************/
class SlowBondset:public FastBondset
{
	//! Each Bondvector is represented by a std::map <VectorInt3, int32_t>, where the identifier is stored as int32_t
	using FastBondset::BondVectors;

	//! True if bondsetLookup is already synchronized with added allowed bond-vectors, false otherwise.
	using FastBondset::lookupSynchronized;

public:

	//! Default constructor
	SlowBondset();

	//! Copy constructor
	SlowBondset(const SlowBondset& copyBondSet);

	//! Default destructor
	virtual ~SlowBondset();

	//! Adding bond to the allowed set
	void addBond(VectorInt3 bondVector, int32_t identifier);

	//! Adding bond to the allowed set
	void addBond(int32_t x, int32_t y, int32_t z, int32_t identifier);

	//! Updates the look-up table of bond-vectors
	void updateLookupTable();

	//! Resets/delete the look-up table of bond-vectors
	void resetLookupTable();

	//! Check if a vector is a valid bond-vector (i.e. part of the set)
	bool isValid(const VectorInt3&) const;

	//! Check if a vector is a valid bond-vector (i.e. part of the set)
	bool isValidStrongCheck(const VectorInt3& bondVector ) const;

private:

  //! lookup table telling if a certain bondvector is valid
  bool*** bondsetLookup;

  //! Coordinates in the bondsetLookup are shifted by this offset
  int32_t lookupOffset;


};

#endif
