/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2016 by
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

#ifndef LEMONADE_UPDATER_MOVES_MOVEADDMONOMERBCC_H
#define LEMONADE_UPDATER_MOVES_MOVEADDMONOMERBCC_H

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBase.h>

/**
 * @file
 * @date   2016/03/14
 * @author Martin
 * @class MoveAddMonomerBcc
 *
 * @brief Add monomer move for the bccBFM to add a vertex/monomer
 *
 **/
template<class TagType=int32_t>
class MoveAddMonomerBcc:public MoveAddMonomerBase<MoveAddMonomerBcc<TagType>, TagType>
{
public:
  MoveAddMonomerBcc(){};
  virtual ~MoveAddMonomerBcc(){};

  //! Reset the probability
  template <class IngredientsType> void init(const IngredientsType& ing);

  //! Check if the move is allowed by the system given as argument.
  template <class IngredientsType> bool check(IngredientsType& ing);

  //! Apply the move to the system given as argument
  template< class IngredientsType> void apply(IngredientsType& ing);
};



/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/**
 * @brief reset the probability and calculate particle index
 */
template <class TagType>
template <class IngredientsType>
void MoveAddMonomerBcc<TagType>::init(const IngredientsType& ing)
{
    this->resetProbability();
    this->setMonomerIndex(ing.getMolecules().size());
}

/*****************************************************************************/
/**
 * @brief check if the move is allowed by the system given as argument.
 */
template <class TagType>
template <class IngredientsType>
bool MoveAddMonomerBcc<TagType>::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}

/*****************************************************************************/
/**
 * @brief apply the move to the system given as argument
 */
template <class TagType>
template <class IngredientsType>
void MoveAddMonomerBcc<TagType>::apply(IngredientsType& ing)
{
  //first add the new monomer at the desired position. this is because
  //some features may want to do things with it
  ing.modifyMolecules().addMonomer(this->getPosition().getX(),this->getPosition().getY(),this->getPosition().getZ());
  this->setMonomerIndex(ing.getMolecules().size()-1);
  //now apply it to the features so that the features can make alterations,
  //for example set the attribute tag, if the FeatureAttributes is used
  ing.applyMove(ing,*this);
}

#endif //LEMONADE_UPDATER_MOVES_MOVEADDMONOMERBCC_H
