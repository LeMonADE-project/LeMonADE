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

#ifndef LEMONADE_FEATURE_FEATUREBONDSETUNSAVECHECK_H
#define LEMONADE_FEATURE_FEATUREBONDSETUNSAVECHECK_H

#include <iostream>
#include <sstream>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/utility/FastBondset.h>
#include <LeMonADE/utility/SlowBondset.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalBase.h>
#include <LeMonADE/updater/moves/MoveConnectBase.h>

/*****************************************************************/
/**
 * @brief Feature adding a set of bond-vectors (Bondset or SlowBondset) to the system.
 * @details Allows bond vectors between refolded monomers!!
 *
 **/
template< class BondSetType=FastBondset>
class FeatureBondsetUnsaveCheck : public FeatureBondset<BondSetType>
{
 public:
	//! Standard constructor (empty)
  FeatureBondsetUnsaveCheck(){}
    
  //! Standard destructor (empty)
  virtual ~FeatureBondsetUnsaveCheck(){}

  using FeatureBondset<BondSetType>::bondset;
  // this is neccessary because we overwrite one checkMove in the following and 
  // therefore we need to explicitly load the other (overloaded ) checkMove-function.
  using FeatureBondset<BondSetType>::checkMove;
  /**
   * @brief Overloaded for MoveConnectBase. 
   *
   * @details Checks if the new bond for this move of type ConnectMoveType is valid.
   * Returns if move is allowed (\a true ) or rejected (\a false ).
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system.
   * @param [in] move A reference to ConnectMoveType.
   * @return if move is allowed (true) or rejected (false).
   */
  template<class IngredientsType,class ConnectMoveType>
  bool checkMove(const IngredientsType& ingredients, const MoveConnectBase<ConnectMoveType>& move) const
  {

	  //get the number of bond partners of the particle to be moved
          uint32_t MonID=move.getIndex();
	  uint32_t partnerID=move.getPartner();
          const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();

	  if (!bondset.isValid(molecules[partnerID]-molecules[MonID])) return false;

          return true;
  }
  
//   using FeatureBondset<BondSetType>::checkMove;
  /**
   * @brief Updates the bond-set lookup table if necessary
   *
   * @details Synchronizes the bond-set look-up table with the rest of the system, and checks for invalid bonds. Connection between periodic images is allowed.
   *
   * @param ingredients A reference to the IngredientsType - mainly the system.
   */
  template<class IngredientsType> void synchronize(IngredientsType& ingredients)
  {
    //this function only does something if the bondset has change since the last update
    bondset.updateLookupTable();
    
    const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
    
    for (size_t i=0; i< molecules.size(); ++i)
    {
     for (size_t j=0; j< molecules.getNumLinks(i); ++j){
       
	 uint n = molecules.getNeighborIdx(i,j);
	 if (!bondset.isValid(molecules[n]-molecules[i]))
	{
	  std::ostringstream errorMessage;
	  errorMessage << "FeatureBondsetUnsaveCheck::synchronize(): Invalid bond vector between monomer " << i << " at " << molecules[i] << " and " << n << " at " <<  molecules[n] <<  ".\n";throw std::runtime_error(errorMessage.str());
	}

    }
    }
  }

 private:

};

#endif
