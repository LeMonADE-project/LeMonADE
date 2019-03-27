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

#ifndef LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUMEBCC_H
#define LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUMEBCC_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>


/*****************************************************************************/
/**
 * @file
 * @date   2016/02/18
 * @author Ron and Martin
 *
 * @class FeatureExcludedVolumeBcc
 * @brief This Feature adds excluded volume check to the system.
 * It is specialised to the body centered cubic lattice.
 *
 * @details The excluded volume is checked by a lattice occupation algorithm
 * on a body centered cubic lattice.
 * The feature is implemented as a class template, where the template parameters
 * \a LatticeClassType and \a LatticeValueTypes specify, what
 * type of lattice is used and which kind of value is saved on the lattice.
 * This feature requires FeatureLattice.
 * Notice: only implementations of Bcc moves included to avoid a mix of
 * FeatureExcludedVolumeBcc and MoveLocalSc.
 *
 * @tparam <LatticeClassType<LatticeValueType>> name of the specialized class.
 * Default is FeatureLattice.
 *
 * @tparam <LatticeValueType> type of the lattice value. Default is bool.
 * */
/*****************************************************************************/

/**
 * @brief Forward declaration. Equivalent to FeatureExcludedVolumeBcc< FeatureLattice <bool> >
 */
template<class LatticeClassType = FeatureLattice<> > class FeatureExcludedVolumeBcc;

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

template<template<typename> class LatticeClassType, typename LatticeValueType>
class FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> > : public Feature {
public:
	//! This Feature requires a lattice.
	typedef LOKI_TYPELIST_1(LatticeClassType<LatticeValueType>) required_features_front;

	//constructor
	FeatureExcludedVolumeBcc() :
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

	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const;

	//! check move for bcc local move
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const;

	//! check sc move: Throw error if wrong lattice Type is used
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalSc& move) const;

	//! check move for adding a bcc monomer
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move) const;

	//! check addsc move: Throw error if wrong lattice Type is used
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move) const;

	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move);

	//! apply move for local bcc moves
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalBcc& move);

	//! apply move for adding an bcc monomer
	template<class IngredientsType, class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerBcc<TagType>& move);

	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

private:

	//! Populates the lattice using the coordinates of molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;

};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveBase& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
//template<class  LatticeClassType<LatticeValueType> >
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove(
		const IngredientsType& ingredients, const MoveBase& move) const
		{
	return true;
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveBase& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing,
const MoveBase& move)
{

}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalSc& move )const
 * @brief Throws a runtime error because the lattice type is inconsitent with the move type
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveLocalSc
 * @return false, throws exception
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove(
const IngredientsType& ingredients, const MoveLocalSc& move) const
{
	throw std::runtime_error("*****FeatureExcludedVolumeSc::check MoveLocalSc: wrong lattice type ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveAddMonomerSc& move )const
 * @brief Throws a runtime error because the lattice type is inconsitent with the move type
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveAddMonomerSc
 * @return false, throws exception
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType, class TagType>
bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove(
const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move) const
{
	throw std::runtime_error("*****FeatureExcludedVolumeSc::check MoveAddMonomerSc: wrong lattice type ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveLocalBcc& move )const
 * @brief checks excluded volume for moves of type MoveLocalBcc
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalBcc.
 * @return if move is allowed (\a true) or rejected (\a false).
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const
{
	if(!latticeFilledUp)
	  throw std::runtime_error("*****FeatureExcludedVolumeBcc::checkMove....lattice is not populated. Run synchronize!\n");

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

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalBcc& move)
 * @brief Updates the lattice occupation according to the move for moves of type MoveLocalBcc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSBcc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalBcc& move)
{
	if (!latticeFilledUp)
		throw std::runtime_error("*****FeatureExcludedVolumeBcc::applyMove....lattice is not populated. Run synchronize!\n");
	//get old position and direction of the move
	VectorInt3 oldPos = ing.getMolecules()[move.getIndex()];
	VectorInt3 direction = move.getDir();

	//change lattice occupation accordingly
	ing.moveOnLattice(oldPos, oldPos + direction);

}

/******************************************************************************/
/**
 * @fn bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveAddMonomerBcc& move )const
 * @brief check excluded volume for insertion of a monomer on a bcc lattice
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerBcc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template < class IngredientsType,class TagType>
bool FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::checkMove( const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move ) const
{
  if (!latticeFilledUp)
	throw std::runtime_error("*****FeatureExcludedVolume_T::checkMove....lattice is not populated. Run synchronize!\n");

  //check if position coordinates are all even or all odd
  VectorInt3 pos=move.getPosition();
  if( ( ((pos.getX())&1) != ((pos.getY())&1) )  || ( ((pos.getX())&1) != ((pos.getZ())&1) ) ){
    return false;
  }

  VectorInt3 dx(1,0,0);
  VectorInt3 dy(0,1,0);
  VectorInt3 dz(0,0,1);

  if(ingredients.getLatticeEntry(pos) ||
    ingredients.getLatticeEntry(pos + dx + dy + dz) || /* check "o" lattice */
    ingredients.getLatticeEntry(pos - dx + dy + dz) ||
    ingredients.getLatticeEntry(pos + dx - dy + dz) ||
    ingredients.getLatticeEntry(pos - dx - dy + dz) ||
    ingredients.getLatticeEntry(pos + dx + dy - dz) ||
    ingredients.getLatticeEntry(pos - dx + dy - dz) ||
    ingredients.getLatticeEntry(pos + dx - dy - dz) ||
    ingredients.getLatticeEntry(pos - dx - dy - dz) ){
    return false;
  }else{
    //check "x" lattice
    int8_t counter(0);
    if(ingredients.getLatticeEntry(pos+ 2*dx)) counter++;
    if(ingredients.getLatticeEntry(pos- 2*dx)) counter++;
    if(ingredients.getLatticeEntry(pos+ 2*dy)) counter++;
    if(ingredients.getLatticeEntry(pos- 2*dy)) counter++;
    if(ingredients.getLatticeEntry(pos+ 2*dz)) counter++;
    if(ingredients.getLatticeEntry(pos- 2*dz)) counter++;

    //! @todo think about counter for "x" lattice check (<2) ... is it too restrictive?
    if(counter < 2)
      return true;
    else
      return false;
  }
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::applyMove( const IngredientsType& ingredients, const MoveAddMonomerBcc& move )const
 * @brief apply excluded volume for MoveAddMonomerBcc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerBcc.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType, class TagType>
void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveAddMonomerBcc<TagType>& move)
{
  ing.setLatticeEntry(move.getPosition(),1);
}

/**
 * Synchronizes the lattice occupation with the rest of the system.
 * The lattice is filled with information about the occupation.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::synchronize(
		IngredientsType& ingredients)
		{
	//note: the lattice entries are set to 0 before by the
	//synchronize function of FeatureLattice
	std::cout << "FeatureExcludedVolumeBcc::synchronizing lattice occupation...";
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

template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureExcludedVolumeBcc< LatticeClassType<LatticeValueType> >::fillLattice(
		IngredientsType& ingredients)
		{
	const typename IngredientsType::molecules_type& molecules =
			ingredients.getMolecules();


	//copy the lattice occupation from the monomer coordinates
	for (size_t n = 0; n < molecules.size(); n++)
	{
		//check if monomer positions are consistent with the bcc model (all even or all odd)
		if( !( ((molecules[n].getX()&1)==1 && (molecules[n].getY()&1)==1 && (molecules[n].getZ()&1)==1) ||
		  ((molecules[n].getX()&1)==0 && (molecules[n].getY()&1)==0 && (molecules[n].getZ()&1)==0) ) ){
		  throw std::runtime_error("FeatureExcludedVolumeBcc<>::fillLattice(): monomer position not consistent with bcc model");
		}

		// iterates over the monomer cube.
		for ( int x = -1 ; x < 2; ++x)
		for ( int y = -1 ; y < 2; ++y)
		for ( int z = -1 ; z < 2; ++z)
		{
			VectorInt3 pos = molecules[n]+VectorInt3(x,y,z);

			if (ingredients.getLatticeEntry(pos))
			{
				std::ostringstream errorMessage;
				errorMessage << "FeatureExcludedVolumeBcc<>::fillLattice(): lattice already occupied when trying to write monomer ";
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
