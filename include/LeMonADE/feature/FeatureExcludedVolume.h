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

#ifndef LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUME_H
#define LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUME_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>


/**
 * @file
 *
 **/

/**
 * @brief Forward declaration. Equivalent to FeatureExcludedVolume< FeatureLattice <bool> >
 */
template<class SpecializedClass = FeatureLattice<> > class FeatureExcludedVolume;

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////


/******************************************************************************/
/**
 *
 *
 *
 *
 * @brief This Feature adds excluded volume check to the system.
 *
 * @details Works with MoveLocalSc and MoveLocalBcc, i.e. it can be used for 
 * both classical bfm and bcc bfm. For the lattice occupation, only a single
 * point is saved on the lattice, in order to be able to use the same feature
 * for both bfm types. This means, that for the excluded volume check 9 lattice 
 * sites need to be checked in classic sc bfm. The feature is implemented as 
 * a class template, where the template parameters \a SpecializedClass and \a ValueTypes specify, what
 * type of lattice is used and which kind of value is saved on the lattice.
 * The default value for \a ValueTypes is bool and the default template for \a SpecializedClass is FeatureLattice
 * If for example some other feature needs to write monomer types on the
 * lattice, one could use FeatureExcludedVolume< FeatureLattice <int8_t> > or
 * for specialized (more performance) a lattice of power 2 type: FeatureExcludedVolume< FeatureLatticePowerOfTwo <int8_t> >.
 * This feature requires FeatureLattice.
 *
 * @tparam <ValueType> type of the lattice value. Default is bool.
 *
 * @tparam <SpecializedClass> name of the specialized class. Default is FeatureLattice.
 **/
template<template<typename> class SpecializedClass, typename ValueType>
class FeatureExcludedVolume< SpecializedClass<ValueType> > : public Feature {
public:
	//! This Feature requires a lattice.
	typedef LOKI_TYPELIST_1(SpecializedClass<ValueType>) required_features_front;

	//constructor
	FeatureExcludedVolume() :
			latticeFilledUp(false) 
	{
	}

	/**
	 * Returns true if the underlying lattice is synchronized and all excluded volume condition
	 * (e.g. monomer/vertex occupies lattice edges) is applied.
	 * Returns false if this feature is out-of-sync.
	 *
	 *
	 * @return true if this feature is synchronized
	 * 		   false if this feature is out-of-sync.
	 **/
	bool isLatticeFilledUp() const {
		return latticeFilledUp;
	}

	/**
	 * Set's the need of synchronization of this feature e.g. escp. if the underlying lattice needs
	 * to refilled and if all excluded volume condition needs to be updated.
	 *
	 *
	 * @param[in] latticeFilledUp Specified if ExVol should be refilled (false) or everything is in-sync (true).
	 *
	 **/
	void setLatticeFilledUp(bool latticeFilledUp) {
		this->latticeFilledUp = latticeFilledUp;
	}

	//! For all unknown moves: this does nothing
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const;

	//! Overloaded for MoveLocalSc (classic ScBFM moves)
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients,
			const MoveLocalSc& move) const;

	//! Overloaded for MoveLocalBcc (BccBFM moves)
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients,
			const MoveLocalBcc& move) const;

	//! For all unknown moves: this does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move);

	//! Overloaded for moves of type MoveLocalBase (applies to both classical ScBFM and BccBFM)
	template<class IngredientsType, class LocalMoveType>
	void applyMove(IngredientsType& ing,
			const MoveLocalBase<LocalMoveType>& move);

	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

private:

	//! Populates the lattice using the (relative) coordinates of the molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;

};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////


/**
 * Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
//template<class LatticeType>
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
bool FeatureExcludedVolume< SpecializedClass<ValueType> >::checkMove(
		const IngredientsType& ingredients, const MoveBase& move) const
		{
	return true;
}


/**
 * Checks excluded volume for moves of type MoveLocalSc. Returns if move is allowed (\a true ) or rejected (\a false ).
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * @return if move is allowed (true) or rejected (false).
 */
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
bool FeatureExcludedVolume< SpecializedClass<ValueType> >::checkMove(
		const IngredientsType& ingredients, const MoveLocalSc& move) const
		{
	if (!latticeFilledUp)
		throw std::runtime_error(
				"*****FeatureExcludedVolume_T::checkMove....lattice is not populated. Run synchronize!\n");

	//get the position of the monomer to be moved (assume "lower left corner")
	VectorInt3 refPos = ingredients.getMolecules()[move.getIndex()];
	//get the direction of the move
	VectorInt3 direction = move.getDir();

	//shift refPos in the direction of the move. this defines then the 
	//point around which the volume occupation must be checked
	refPos += (2 * direction);

	//now the nine closest positions in the plane of refPos perpendicular to the 
	//move direction must be checked

	/*get two directions perpendicular to vector direction of the move*/
	VectorInt3 perp1, perp2;
	/* first perpendicular direction is either (0 1 0) or (1 0 0)*/
	int32_t x1 = ((direction.getX() == 0) ? 1 : 0);
	int32_t y1 = ((direction.getX() != 0) ? 1 : 0);
	perp1.setX(x1);
	perp1.setY(y1);
	perp1.setZ(0);

	/* second perpendicular direction is either (0 0 1) or (0 1 0)*/
	int32_t y2 = ((direction.getZ() == 0) ? 0 : 1);
	int32_t z2 = ((direction.getZ() != 0) ? 0 : 1);
	perp2.setX(0);
	perp2.setY(y2);
	perp2.setZ(z2);

	VectorInt3 v1 = perp1 + perp2;
	VectorInt3 v2 = perp1 - perp2;
	//check if the lattice sites are free
	if (
	(ingredients.getLatticeEntry(refPos)) ||
			(ingredients.getLatticeEntry(refPos + perp1)) ||
			(ingredients.getLatticeEntry(refPos - perp1)) ||
			(ingredients.getLatticeEntry(refPos + perp2)) ||
			(ingredients.getLatticeEntry(refPos - perp2)) ||
			(ingredients.getLatticeEntry(refPos + v1)) ||
			(ingredients.getLatticeEntry(refPos - v1)) ||
			(ingredients.getLatticeEntry(refPos - v2)) ||
			(ingredients.getLatticeEntry(refPos + v2))
			)
		return false;
	else
		return true;

}


/**
 * Checks excluded volume for moves of type MoveLocalBcc. Returns if move is allowed (\a true ) or rejected (\a false ).
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalBcc.
 * @return if move is allowed (true) or rejected (false).
 */
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
bool FeatureExcludedVolume< SpecializedClass<ValueType> >::checkMove(
		const IngredientsType& ingredients, const MoveLocalBcc& move) const
		{
#ifdef DEBUG
	if(!latticeFilledUp)
	throw std::runtime_error("*****FeatureExcludedVolume_T::checkMove....lattice is not populated. Run synchronize!\n");
#endif

	int32_t x, y, z;
	int8_t dx, dy, dz;
	x = ingredients.getMolecules()[move.getIndex()][0];
	y = ingredients.getMolecules()[move.getIndex()][1];
	z = ingredients.getMolecules()[move.getIndex()][2];

	dx = 2 * move.getDir()[0];
	dy = 2 * move.getDir()[1];
	dz = 2 * move.getDir()[2];

	if (ingredients.getLatticeEntry(x + dx, y, z) ||
			ingredients.getLatticeEntry(x, y + dy, z) ||
			ingredients.getLatticeEntry(x, y, z + dz) ||
			ingredients.getLatticeEntry(x + dx, y + dy, z) ||
			ingredients.getLatticeEntry(x + dx, y, z + dz) ||
			ingredients.getLatticeEntry(x, y + dy, z + dz) ||
			ingredients.getLatticeEntry(x + dx, y + dy, z + dz)
					)
		return false;
	else
		return true;
}

/**
 * This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
void FeatureExcludedVolume< SpecializedClass<ValueType> >::applyMove(IngredientsType& ing,
		const MoveBase& move)
		{

}


/**
 * This function updates the lattice occupation according to the move (for moves of type MoveLocalBase ).
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType, class LocalMoveType>
void FeatureExcludedVolume< SpecializedClass<ValueType> >::applyMove(IngredientsType& ing,
		const MoveLocalBase<LocalMoveType>& move)
		{
	if (!latticeFilledUp)
		throw std::runtime_error(
				"*****FeatureExcludedVolume_T::applyMove....lattice is not populated. Run synchronize!\n");
	//get old position and direction of the move
	VectorInt3 oldPos = ing.getMolecules()[move.getIndex()];
	VectorInt3 direction = move.getDir();

	//change lattice occupation accordingly
	ing.moveOnLattice(oldPos, oldPos + direction);

}


/**
 * Synchronizes the lattice occupation with the rest of the system.
 * The lattice is filled with information about the occupation.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
void FeatureExcludedVolume< SpecializedClass<ValueType> >::synchronize(
		IngredientsType& ingredients)
		{
	//note: the lattice entries are set to 0 before by the 
	//synchronize function of FeatureLattice
	std::cout << "FeatureExcludedVolume_T::synchronizing lattice occupation...";
	fillLattice(ingredients);
	std::cout << "done\n";
}


/**
 * This function populates the lattice directly with positions from molecules
 * e.g. sets the value at the (relative) monomer position on the lattice to 1.
 * It also has a simple check if the target lattice is already occupied.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */

template<template<typename> class SpecializedClass, typename ValueType>
template<class IngredientsType>
void FeatureExcludedVolume< SpecializedClass<ValueType> >::fillLattice(
		IngredientsType& ingredients)
		{
	const typename IngredientsType::molecules_type& molecules =
			ingredients.getMolecules();
			

	//copy the lattice occupation from the monomer coordinates
	for (size_t n = 0; n < molecules.size(); n++)
	{
		// iterates over the monomer cube.

		for ( int x = -1 ; x < 2; ++x)
		for ( int y = -1 ; y < 2; ++y)
		for ( int z = -1 ; z < 2; ++z)
		{		      
			VectorInt3 pos = molecules[n]+VectorInt3(x,y,z);
			
			if (ingredients.getLatticeEntry(pos))
			{
				std::ostringstream errorMessage; 
				errorMessage << "FeatureExcludedVolume<>::fillLattice(): lattice already occupied when trying to write monomer ";
				errorMessage << n << " at " << molecules[n] << ", cube corner at " << pos << ".\n"; 
				throw std::runtime_error(errorMessage.str());
			}
			//here we simply set a one on every occupied lattice
			//site. this assumes that 1 can be cast to LatticeType,
			//even though in principle LatticeType could be anything.
			//note that this may be just a preliminiary initialization,
			//as other features may assign more specific values to the
			//lattice site.
		}
		ingredients.setLatticeEntry(molecules[n],1);
	}
	latticeFilledUp = true;
}


#endif
