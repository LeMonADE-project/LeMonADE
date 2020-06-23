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

#ifndef LEMONADE_FEATURE_FEATUREBREAK_H
#define LEMONADE_FEATURE_FEATUREBREAK_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalBcc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBcc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>
#include <LeMonADE/utility/Lattice.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
/*****************************************************************************/
/**
 * @file
 * @date   2020/06/06
 * @author Toni
 *
 * @class FeatureBreak
 * @brief This Feature erases bonds between monomers.
 *
 * @details Works for all lattice types. Pretty simple implementation that 
 * only checks if there is a bond betwenn 
 *
 * @tparam 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureBreak : public Feature {
  
public:
	//! check move for basic move - always true
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const{ return true;};
	
    //! check bas connect move - always true 
	template<class IngredientsType, class SpecializedMove> 
	bool checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const;

    //! check bas connect move - always true 
	template<class IngredientsType > 
	bool checkMove(const IngredientsType& ingredients, const MoveBreak& move) const;
};

///////////////////////////////////////////////////////////////////////////////
////////////////////////// member definitions /////////////////////////////////
/******************************************************************************/
/**
 * @fn void FeatureBreak ::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const
 * @brief This function checks if the monomers are connected at all. 
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveBreakBase is a base class for all breaking moves.
 */
/******************************************************************************/
template<class IngredientsType, class SpecializedMove> 
bool FeatureBreak::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const
{
  if ( ! ingredients.getMolecules().areConnected(move.getIndex(),move.getPartner())) return false; 
  return true;
}
/******************************************************************************/
/**
 * @fn void FeatureBreak ::checkMove(const IngredientsType& ingredients, const MoveBreakBase<SpecializedMove>& move) const
 * @brief This functions returns always true, because the move itself serves IDs which are connected.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system
 * @param [in] move MoveBreak breaks a bond in the system. 
 */
/******************************************************************************/
template<class IngredientsType > 
bool FeatureBreak::checkMove(const IngredientsType& ingredients, const MoveBreak& move) const
{
  return true;
}
#endif
