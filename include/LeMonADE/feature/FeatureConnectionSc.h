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
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveBreakBase.h>
#include <LeMonADE/updater/moves/MoveBreak.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/utility/Lattice.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>

/*****************************************************************************/
/**
 * @file
 * @date   2019/02/05
 * @author Toni
 *
 * @class FeatureConnectionSc
 * @brief This Feature add new bonds between monomers.
 *
 * @details Works only in combination with an excluded volume feature. 
 * Updated: 2020/06/17
 * The feature creates a lattice where the monomer ID of all monomers are written on. 
 * This is used to find bond partners by the ConnectionMoves. The connection moves 
 * only checks if the maxNumber of links is reached and of the monomers are already
 * connected. The spatial moves update the IDs on the lattice. 
 *
 * @tparam 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureConnectionSc : public Feature {
  
public:
	
	//! this feature will not use any of FeatureExcludedVolumeSc<> but we need excluded volume property
    // 	typedef LOKI_TYPELIST_1(FeatureExcludedVolumeSc<>) required_features_back;

	//constructor
	FeatureConnectionSc() :latticeFilledUp(false)
	{connectionLattice.setupLattice();}
	
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
	
	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const{ return true;};

	//! check bas connect move - always true 
	template<class IngredientsType, class SpecializedMove> 
	bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const;
	
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const {throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalBcc: wrong lattice type ... \n"); return false;};
	
	//! check bas connect move - always true 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move) const {throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalBcc: wrong lattice type ... \n"); return false;};
	
	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move){};
	
	//!apply move for the scBFM local move which changes the lattice
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalSc& move);	
	
	//!apply move for the scBFM local diagonal move which changes the lattice
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalScDiag& move);	
        
	//!
	template<class IngredientsType, class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move);	
	
	//! Synchronize with system: Fill the lattice with 1 (occupied) and 0 (free).
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients);
	
	//! Get the lattice value at a certain point
	//! return monomer index between (0; molecules.size()-1) if there's a monomer
	//! return uint32_t(-1)=4294967295 if place is empty
	uint32_t getIdFromLattice(const VectorInt3& pos) const { return connectionLattice.getLatticeEntry(pos)-1;} ;

	//! Get the lattice value at a certain point
	//! return monomer index between (0; molecules.size()-1) if there's a monomer
	//! return uint32_t(-1)=4294967295 if place is empty
	uint32_t getIdFromLattice(const int x, const int y, const int z) const { return connectionLattice.getLatticeEntry(x,y,z)-1;};

protected:

	//! Populates the lattice using the coordinates of molecules.
	template<class IngredientsType> void fillLattice(
			IngredientsType& ingredients);

	//! Tag for indication if the lattice is populated.
	bool latticeFilledUp;
	//!
	Lattice<uint32_t> connectionLattice;
    
    uint32_t maxConnectivity;
};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc::checkMove( const IngredientsType& ingredients, const MoveConnectSc& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 * @details  it might make a difference for the speed if the order of statements is switched for different systems parameters
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
bool FeatureConnectionSc ::checkMove(const IngredientsType& ingredients, const MoveConnectBase<SpecializedMove>& move) const
{
  
  	if (!latticeFilledUp)
	    throw std::runtime_error("*****FeatureConnectionSc::checkMove....lattice is not populated. Run synchronize!\n");

    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    
    uint32_t ID(move.getIndex());
	//check for maximum number of bonds for the first monomer
	if ( molecules.getNumLinks(ID) >=  maxConnectivity) return false;

	uint32_t Neighbor(move.getPartner());
	//check for maximum number of bonds for the second monomer
	if ( molecules.getNumLinks(Neighbor) >= maxConnectivity ) return false;

	//check if the two monomers are already connected
	if ( molecules.areConnected(ID,Neighbor) ) return false;

	//if still here, then the two monomers are allowed to connect 
	return true;
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalSc& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveLocalSc& move)
{
  VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
  VectorInt3 direction=move.getDir();
  VectorInt3 oldPlusDir=oldPos+direction;
  connectionLattice.moveOnLattice(oldPos,oldPlusDir);
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalScDiag& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveLocalScDiag& move)
{
  VectorInt3 oldPos=ing.getMolecules()[move.getIndex()];
  VectorInt3 direction=move.getDir();
  VectorInt3 oldPlusDir=oldPos+direction;
  connectionLattice.moveOnLattice(oldPos,oldPlusDir);
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move)
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
template<class IngredientsType, class TagType>
void FeatureConnectionSc  ::applyMove(IngredientsType& ing,const MoveAddMonomerSc<TagType>& move)
{
  uint32_t MonID(move.getMonomerIndex()); 
  VectorInt3 pos=ing.getMolecules()[MonID];
  connectionLattice.setLatticeEntry(pos,MonID+1 );
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
template<class IngredientsType>
void FeatureConnectionSc  ::synchronize(IngredientsType& ingredients)
{

	std::cout << "FeatureConnectionSc::synchronizing lattice occupation...\n";
	fillLattice(ingredients);
    maxConnectivity=ingredients.getMolecules().getMaxConnectivity();
	std::cout << "done\n";
}


/******************************************************************************/
/**
 * @fn void FeatureExcludedVolumeSc ::fillLattice(IngredientsType& ingredients)
 * @brief This function populates the lattice directly with positions from molecules.
 * It also has a simple check if the target lattice is already occupied.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 * */
/******************************************************************************/
template<class IngredientsType>
void FeatureConnectionSc::fillLattice(IngredientsType& ingredients)
{
  
	connectionLattice.setupLattice(ingredients.getBoxX(),ingredients.getBoxY(),ingredients.getBoxZ());
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	//copy the lattice occupation from the monomer coordinates
	for(size_t n=0;n<molecules.size();n++)
	{
		VectorInt3 pos=ingredients.getMolecules()[n];
		if( connectionLattice.getLatticeEntry(pos)!=0 )
		{
			throw std::runtime_error("********** FeatureConnectionSc::fillLattice: multiple lattice occupation ******************");
		}
		else 
		{
			//here we simply set the monomer id (plus one!) on the lattice site 
			// the offset implies that the index zero is still used for unoccupied
			// with and unreactive monomer
			VectorInt3 pos=molecules[n];
			connectionLattice.setLatticeEntry(pos,n+1);
		}
	}
	latticeFilledUp=true;
}

#endif
