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

#ifndef LEMONADE_FEATURE_FEATURECONNECTIONSC_H
#define LEMONADE_FEATURE_FEATURECONNECTIONSC_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>

/**
 * @class MonomerReactivity
 * @brief set the monomer reactive or unreactive 
 * 
 */
class MonomerReactivity
{
  public:

	//! Standard constructor- initially the reactivity is set to false.
	MonomerReactivity():reactivity(false),nMaxBonds(2){}

	//! Getting the reactivity of the monomer.
	bool IsReactive() const {return reactivity;}
	
	//! Getting the number of maximum possible bonds for the monomer.
	uint32_t getNMaxBonds() const {return nMaxBonds;};

	/**
	 * @brief Setting the reactivity of the monomer with \para reactivity_.
	 *
	 * @param reactivity_ either trueor false
	 */
	void setReactive(bool reactivity_){ reactivity=reactivity_;}
	/**
	 * @brief Setting the maximum possible bonds of the monomer with \para nMaxBonds_.
	 *
	 * @param nMaxBonds_ 
	 */	
	void setNMaxBonds(uint32_t nMaxBonds_){nMaxBonds=nMaxBonds_;} 

private:
     //! Private variable holding the tag. Default is NULL.
     bool reactivity;
     //!
     uint32_t nMaxBonds;
     
};

/*****************************************************************************/
/**
 * @file
 * @date   2019/02/05
 * @author Toni
 *
 * @class FeatureConnectionSc
 * @brief This Feature add new bonds between monomers.
 *
 * @details 
 *
 * @tparam 
 * */
/*****************************************************************************/
/**
 * @brief Forward declaration. Equivalent to FeatureExcludedVolumeSc< FeatureLattice <uint32_t> >
 */
template<class LatticeClassType = FeatureLattice<uint32_t> > class FeatureConnectionSc;

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////
template<template<typename> class LatticeClassType, typename LatticeValueType>
class FeatureConnectionSc< LatticeClassType<LatticeValueType> > : public Feature {
  
public:
  	//! This Feature requires a lattice.
	typedef LOKI_TYPELIST_1(LatticeClassType<LatticeValueType>) required_features_front;
	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerReactivity) monomer_extensions;

	/**
	 * Returns true if the underlying lattice is synchronized and all excluded volume condition
	 * (e.g. monomer/vertex occupies lattice edges) is applied.
	 * Returns false if this feature is out-of-sync.
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
	 * @param[in] latticeFilledUp Specified if ExVol should be refilled (false) or everything is in-sync (true).
	 *
	 **/
	void setLatticeFilledUp(bool latticeFilledUp) {
		this->latticeFilledUp = latticeFilledUp;
	}

	//constructor
	FeatureConnectionSc() :latticeFilledUp(false)
	{}

// 	//! check move for basic move - always true
// 	template<class IngredientsType>
// 	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const;

	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveConnectBase& move) const;
	
	//! check bas connect move 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveConnectSc& move) const;
	
	
// 	//! apply move for basic moves - does nothing
// 	template<class IngredientsType>
// 	void applyMove(IngredientsType& ing, const MoveBase& move);
	
	//! apply move for the baseic connection - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveConnectBase& move) const;
	
	//!apply move for the scBFM connection move for connection
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveConnectSc& move) const;
	
	
	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);

protected:

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
 * @fn bool FeatureConnectionSc::checkMove( const IngredientsType& ingredients, const MoveConnectBase& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureConnectionSc::checkMove(const IngredientsType& ingredients, const MoveConnectBase& move) const
{
	return true;
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveConnectBase& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureConnectionSc< LatticeClassType<LatticeValueType> > ::applyMove(IngredientsType& ing,const MoveConnectBase& move)
{

}
/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc::checkMove( const IngredientsType& ingredients, const MoveConnectSc& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
bool FeatureConnectionSc< LatticeClassType<LatticeValueType> >::checkMove(const IngredientsType& ingredients, const MoveConnectSc& move) const
{
  ///@todo implementation 
  
	return true;
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveConnectSc& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureConnectionSc< LatticeClassType<LatticeValueType> > ::applyMove(IngredientsType& ing,const MoveConnectSc& move)
{

}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the lattice occupation with the rest of the system
 * by calling the private function fillLattice.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureConnectionSc< LatticeClassType<LatticeValueType> > ::synchronize(IngredientsType& ingredients)
{

	std::cout << "FeatureConnectionSc::synchronizing lattice occupation...\n";
	fillLattice(ingredients);
	std::cout << "done\n";
}


/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc< LatticeClassType<LatticeValueType> >::fillLattice(IngredientsType& ingredients)
 * @brief This function populates the lattice directly with positions from molecules.
 * It also has a simple check if the target lattice is already occupied.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 * */
/******************************************************************************/
template<template<typename> class LatticeClassType, typename LatticeValueType>
template<class IngredientsType>
void FeatureConnectionSc< LatticeClassType<LatticeValueType> >::fillLattice(IngredientsType& ingredients)
{
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	//copy the lattice occupation from the monomer coordinates
	for(size_t n=0;n<molecules.size();n++)
	{
		VectorInt3 pos=ingredients.getMolecules()[n];
		if( ingredients.getLatticeEntry(pos)!=0 )
		{
			throw std::runtime_error("********** FeatureConnectionSc::fillLattice: multiple lattice occupation ******************");
		}
		else if (ingredients.getMolecules()[n].IsReactive())
		{
			//here we simply set the monomer id (plus one!) on the lattice site 
			// the offset implies that the index zero is still used for unoccupied
			// with and unreactive monomer
			VectorInt3 pos=molecules[n];
			ingredients.setLatticeEntry(pos,n+1);
		}
	}
	latticeFilledUp=true;
}
#endif
