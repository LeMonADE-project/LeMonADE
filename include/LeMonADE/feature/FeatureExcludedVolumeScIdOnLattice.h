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

#ifndef LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUMESC_IDONLATTICE_H
#define LEMONADE_FEATURE_FEATUREEXCLUDEDVOLUMESC_IDONLATTICE_H

#include <LeMonADE/feature/FeatureLattice.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/utility/LatticePredicates.h>

/*****************************************************************************/
/**
 * @file
 * @date   2017/09/19
 * @author Toni
 *
 * @class FeatureExcludedVolumeScIdOnLattice
 * @brief This Feature adds excluded volume check to the system.
 * It is specialised to the simple cubic lattice and write a value using the predicate 
 * on lattice.
 *
 * @details The excluded volume is checked by a lattice occupation algorithm
 * on a simple cubic lattice.
 * The feature is implemented as a class template, where the template parameters
 * \a LatticeClassType and \a LatticeValueTypes specify, what
 * type of lattice is used and which kind of value is saved on the lattice.
 * This feature requires FeatureLattice.
 * Notice: only implementations of Sc moves included to avoid a mix of
 * FeatureExcludedVolumeScIdOnLattice and MoveLocalBcc. 
 * To use the feature within the LOKI_TYPELIST_NN first declare something like :
 * typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, MonomerID  > MyExcludedVolumeFeature;
 * and use the MyExcludedVolumeFeature for the typelist. Otherwise it appears as two arguments for the Loki.
 *
 * @tparam <LatticeClassType<LatticeValueType>> name of the specialized class.
 * Default is FeatureLattice.
 *
 * @tparam <LatticeValueType> type of the lattice value. Default is bool.
 * */
/*****************************************************************************/

/**
 * @brief Forward declaration. Equivalent to FeatureExcludedVolumeScIdOnLattice< FeatureLattice <bool> >
 */

template<class LatticeClassType = FeatureLatticePowerOfTwo<uint32_t>, class Predicate = Bool  > class FeatureExcludedVolumeScIdOnLattice;

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

template<template<typename> class LatticeClassType, typename LatticeValueType, class Predicate>
class FeatureExcludedVolumeScIdOnLattice<LatticeClassType<LatticeValueType>, Predicate > : public  FeatureExcludedVolumeSc<LatticeClassType<LatticeValueType> > {

	typedef FeatureExcludedVolumeSc<LatticeClassType<LatticeValueType> > BaseClass;
	using BaseClass::latticeFilledUp;

public:	
	//This feature uses the attributes to write on the lattice for some predicates
	typedef LOKI_TYPELIST_1(FeatureAttributes<>) required_features_back;

	//constructor
	FeatureExcludedVolumeScIdOnLattice(const Predicate& pred_ = Predicate()) :
			BaseClass(),
			pred(pred_)
	{
	}
	
	using BaseClass::applyMove;

	//! apply move for adding an sc monomer
	template<class IngredientsType,class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move);
	
	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
	

protected:
	//! Populates the lattice using the coordinates of molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);
private:
  
	//! predicate for the lattice occupation 
	Predicate pred; 
};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/******************************************************************************/

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeScIdOnLattice< LatticeClassType<LatticeValueType> >::applyMove( const IngredientsType& ingredients, const MoveAddMonomerSc& move )const
 * @brief apply excluded volume for MoveAddMonomerSc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerSc.
 * */
/******************************************************************************/
// template<template<typename> class LatticeClassType, typename LatticeValueType>
template<template<typename> class LatticeClassType, typename LatticeValueType, class Predicate>
template<class IngredientsType, class TagType>
void FeatureExcludedVolumeScIdOnLattice< LatticeClassType<LatticeValueType> , Predicate>::applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move)
{
  VectorInt3 pos=move.getPosition();
  VectorInt3 dx(1,0,0);
  VectorInt3 dy(0,1,0);
  VectorInt3 dz(0,0,1);
 std::cout <<"use the correct appy move "<<std::endl;
  LatticeValueType value=move.getMonomerIndex()+1;
  ing.setLatticeEntry(pos,value);
  ing.setLatticeEntry(pos+dx,value);
  ing.setLatticeEntry(pos+dy,value);
  ing.setLatticeEntry(pos+dx+dy,value);
  ing.setLatticeEntry(pos+dz,value);
  ing.setLatticeEntry(pos+dz+dx,value);
  ing.setLatticeEntry(pos+dz+dy,value);
  ing.setLatticeEntry(pos+dz+dx+dy,value);
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeScIdOnLattice< LatticeClassType<LatticeValueType> >::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the lattice occupation with the rest of the system
 * by calling the private function fillLattice.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
// template<template<typename> class LatticeClassType, typename LatticeValueType>
template<template<typename> class LatticeClassType, typename LatticeValueType, class Predicate>
template<class IngredientsType>
void FeatureExcludedVolumeScIdOnLattice< LatticeClassType<LatticeValueType>, Predicate >::synchronize(IngredientsType& ingredients)
{
	//note: the lattice entries are set to 0 before by the
	//synchronize function of FeatureLattice
	std::cout << "FeatureExcludedVolumeScIdOnLattice::synchronizing lattice occupation...\n";
	fillLattice(ingredients);
	std::cout << "done\n";
}

/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeScIdOnLattice< LatticeClassType<LatticeValueType> >::fillLattice(IngredientsType& ingredients)
 * @brief This function populates the lattice directly with positions from molecules.
 * It also has a simple check if the target lattice is already occupied.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 * */
/******************************************************************************/
// template<template<typename> class LatticeClassType, typename LatticeValueType>
template<template<typename> class LatticeClassType, typename LatticeValueType, class Predicate>
template<class IngredientsType>
void FeatureExcludedVolumeScIdOnLattice< LatticeClassType<LatticeValueType>, Predicate >::fillLattice(IngredientsType& ingredients)
{
  
  std::cout << "Use the correct fillLattice" << std::endl;
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	//copy the lattice occupation from the monomer coordinates
	for(size_t n=0;n<molecules.size();n++)
	{
		VectorInt3 pos=ingredients.getMolecules()[n];

		if( ingredients.getLatticeEntry(pos)!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,0,0))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(0,1,0))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(0,0,1))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,1,0))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,0,1))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(0,1,1))!=0 ||
		  ingredients.getLatticeEntry(pos+VectorInt3(1,1,1))!=0
		)
		{
			throw std::runtime_error("********** FeatureExcludedVolume::fillLattice: multiple lattice occupation ******************");
		}
		else
		{
			//here we use the functor pred to decide what should be written on the lattice
			VectorInt3 pos=molecules[n];
			LatticeValueType value=pred(molecules,n);
			ingredients.setLatticeEntry(pos,value);
			ingredients.setLatticeEntry(pos+VectorInt3(1,0,0),value);
			ingredients.setLatticeEntry(pos+VectorInt3(0,1,0),value);
			ingredients.setLatticeEntry(pos+VectorInt3(1,1,0),value);
			ingredients.setLatticeEntry(pos+VectorInt3(0,0,1),value);
			ingredients.setLatticeEntry(pos+VectorInt3(1,0,1),value);
			ingredients.setLatticeEntry(pos+VectorInt3(0,1,1),value);
			ingredients.setLatticeEntry(pos+VectorInt3(1,1,1),value);
		}

	}
	latticeFilledUp=true;
}

#endif
