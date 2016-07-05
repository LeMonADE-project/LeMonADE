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

#ifndef LEMONADE_UPDATER_MOVES_MOVEADDSCMONOMER_H
#define LEMONADE_UPDATER_MOVES_MOVEADDSCMONOMER_H

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveBase.h>

/**
 * @file
 * 
 * @author Hauke
 * @class MoveAddScMonomer
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM to add a vertex/monomer
 *
 * @deprecated Untested!!
 *
 * @todo Test!
 * 
 * @todo Code doubling with MoveAddBccMonomer! Specialized class only neccessary to enable the features to handle the moves types differently.
 **/
class MoveAddScMonomer:public MoveBase
{
public:
  
  //! Reset the probability
  template <class IngredientsType> void init(const IngredientsType& ing);

  //! Check if the move is allowed by the system given as argument.
  template <class IngredientsType> bool check(IngredientsType& ing);

  //! Apply the move to the system given as argument
  template< class IngredientsType> void apply(IngredientsType& ing);
  
  //! setter function for type of new monomer
  void setType(int32_t t){type=t;}
  //! setter function for position of new monomer taking a VectorInt3
  void setPosition(VectorInt3 pos){position=pos;}
  //! setter function for position of new monomer taking a triple of ints
  void setPosition(int32_t x,int32_t y,int32_t z){position.setX(x);position.setY(y);position.setZ(z);}
  //! getter function for the type of the new monomer
  int32_t getType() const{return type;}
  //! getter function for the position of the new monomer returning a VectorInt3
  const VectorInt3 getPosition() const {return position;}
  //! getter function for index of the new monomer. This is ing.getMolecules().size() before applying and ing.getMolecules().size()-1 after applying the move.
  size_t getParticleIndex() const {return particleIndex;}
  
private:
  //! position where the new monomer is placed in the simulation box
  VectorInt3 position;
  //! type that is applied to the new monomer, requires Feature FeatureAttributes with int-Type
  int32_t type;
  /** 
   * @brief Index of new PartileThis is ing.getMolecules().size() before applying and ing.getMolecules().size()-1 after applying the move.
   * @details It is set when apply is called. useful if Features want to alter the particle when applying the move
   */
  size_t particleIndex;
};



/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/**
 * @brief reset the probability
 */
template <class IngredientsType>
void MoveAddScMonomer::init(const IngredientsType& ing)
{
    this->resetProbability();
    particleIndex=ing.getMolecules().size();
}

/*****************************************************************************/
/**
 * @brief check if the move is allowed by the system given as argument.
 */
template <class IngredientsType>
bool MoveAddScMonomer::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}
  
/*****************************************************************************/
/**
 * @brief apply the move to the system given as argument
 */
template< class IngredientsType>
void MoveAddScMonomer::apply(IngredientsType& ing)
{
  //first add the new monomer at the desired position. this is because
  //some features may want to do things with it
  ing.modifyMolecules().addMonomer(position.getX(),position.getY(),position.getZ());
  particleIndex=ing.getMolecules().size()-1;
  //now apply it to the features so that the features can make alterations,
  //for example set the attribute tag, if the FeatureAttributes is used
  ing.applyMove(ing,*this);
}

#endif //LEMONADE_UPDATER_MOVES_MOVEADDSCMONOMER_H
