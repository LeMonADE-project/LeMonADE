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
#include <LeMonADE/feature/FeatureConnectionScIdOnLattice.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>

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


///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureConnectionSc : public Feature {
  
typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, MonomerID  > ExcludedVolumeFeature;
public:
  
	//! This Feature requires a lattice.
	typedef LOKI_TYPELIST_1(ExcludedVolumeFeature) required_features_front;

	//constructor
	FeatureConnectionSc() :
	{}

	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const;

	//! check move for sc local move: Throw error if wrong lattice Type is used
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalSc& move) const;

	//! check move for sc local diagonal move : Throw error if wrong lattice Type is used
	template<class IngredientsType>	
	bool checkMove( const IngredientsType& ingredients, const MoveLocalScDiag& move ) const;
	//! check bcc move: Throw error if wrong lattice Type is used
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const;

	//! check move for adding an sc monomer : Throw error if wrong lattice Type is used
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move) const;

	//! check bcc addmove: Throw error if wrong lattice Type is used
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& move) const;

	//! check bas connect move - always true 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveConnectBase& move) const;
	
	//! check bas connect move 
	template<class IngredientsType, class TagType>
	bool checkMove(const IngredientsType& ingredients, const MoveConnectSc& move) const;
	
	//! apply move for basic moves - does nothing
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveBase& move);

	//! apply move for local sc moves
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalSc& move);
	
	//! apply move for local sc diagonal moves 
	template<class IngredientsType>
	void applyMove(IngredientsType& ing, const MoveLocalScDiag& move);
	
	//! apply move for adding an sc monomer
	template<class IngredientsType, class TagType>
	void applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move);

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
 * @fn bool FeatureConnectionSc::checkMove( const IngredientsType& ingredients, const MoveBase& move )const
 * @brief Returns true for all moves other than the ones that have specialized versions of this function.
 * This dummy function is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 * @return true Always!
 */
/******************************************************************************/
template<class IngredientsType>
bool FeatureConnectionSc::checkMove(const IngredientsType& ingredients, const MoveBase& move) const
{
	return true;
}

/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveBase& move)
 * @brief This function applies for unknown moves other than the ones that have specialized versions of this function.
 * It does nothing and is implemented for generality.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move General move other than MoveLocalSc or MoveLocalBcc.
 */
/******************************************************************************/
 
template<class IngredientsType>
void FeatureConnectionSc ::applyMove(IngredientsType& ing,const MoveBase& move)
{

}

/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveLocalBcc& move )const
 * @brief Throws a runtime error because the lattice type is inconsitent with the move type
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveLocalBcc, which is the wrong one in this case
 * @return false, throws an exception
 */
/******************************************************************************/
 
template<class IngredientsType>
bool FeatureConnectionSc ::checkMove(const IngredientsType& ingredients, const MoveLocalBcc& move) const
{
	throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalBcc: wrong move for connection ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveAddMonomerBcc& addmove )const
 * @brief Throws a runtime error because the lattice type is inconsitent with the move type
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveAddMonomerBcc, which is the wrong one in this case
 * @return false, throws an exception
 */
/*****************************************************************************
*/

 
template<class IngredientsType, class TagType>
bool FeatureConnectionSc ::checkMove(const IngredientsType& ingredients, const MoveAddMonomerBcc<TagType>& addmove) const
{
	throw std::runtime_error("*****FeatureConnectionSc::check MoveAddMonomerBcc: wrong move for connection ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveLocalSc& move )const
 * @brief checks excluded volume for moves of type MoveLocalSc
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * @return if move is allowed (\a true) or rejected (\a false).
 * */
/******************************************************************************/
 
template < class IngredientsType>
bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveLocalSc& move ) const
{
	throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalSc: wrong move for connection ... \n");
	return false;
}

/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveLocalScDiag& move )const
 * @brief checks excluded volume for moves of type MoveLocalScDiag
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * @return if move is allowed (\a true) or rejected (\a false).
 * */
/******************************************************************************/
 
template < class IngredientsType>
bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveLocalScDiag& move ) const
{
	throw std::runtime_error("*****FeatureConnectionSc::check MoveLocalScDiag: wrong move for connection ... \n");
	return false;
}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalSc<LocalMoveType>& move)
 * @brief Updates the lattice occupation according to the move for moves of type MoveLocalSc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * */
/******************************************************************************/
 
template<class IngredientsType>
void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalSc& move)
{
	

}
/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveLocalScDiag& move)
 * @brief updates the lattice ocupation according to the move for types MoveLocalScDiag 
 * 
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveLocalSc.
 * */
/******************************************************************************/
 
template<class IngredientsType>
void FeatureConnectionSc<LatticeClassType<LatticeValueType> >::applyMove(IngredientsType& ing, const MoveLocalScDiag& move)
{       
	
}

/******************************************************************************/
/**
 * @fn bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveAddMonomerSc& move )const
 * @brief check excluded volume for insertion of a monomer on a sc lattice
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerSc.
 * */
/******************************************************************************/
 
template < class IngredientsType, class TagType>
bool FeatureConnectionSc ::checkMove( const IngredientsType& ingredients, const MoveAddMonomerSc<TagType>& move ) const
{
 	throw std::runtime_error("*****FeatureConnectionSc::check MoveAddMonomerSc: wrong move for connection ... \n");
	return false;
}


/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::applyMove( const IngredientsType& ingredients, const MoveAddMonomerSc& move )const
 * @brief apply excluded volume for MoveAddMonomerSc.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move A reference to MoveAddMonomerSc.
 * */
/******************************************************************************/
 
template<class IngredientsType, class TagType>
void FeatureConnectionSc ::applyMove(IngredientsType& ing, const MoveAddMonomerSc<TagType>& move)
{
}

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
 
template<class IngredientsType>
void FeatureConnectionSc ::applyMove(IngredientsType& ing,const MoveConnectBase& move)
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
template<class IngredientsType>
bool FeatureConnectionSc::checkMove(const IngredientsType& ingredients, const MoveConnectSc& move) const
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
 
template<class IngredientsType>
void FeatureConnectionSc ::applyMove(IngredientsType& ing,const MoveConnectSc& move)
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
 
template<class IngredientsType>
void FeatureConnectionSc ::synchronize(IngredientsType& ingredients)
{

// 	std::cout << "FeatureConnectionSc::synchronizing lattice occupation...\n";
// 	std::cout << "done\n";
}
#endif
