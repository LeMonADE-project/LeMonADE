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

#include <LeMonADE/updater/moves/MoveBase.h>

/**
 * @class MoveAddScMonomer
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM to add a vertex/monomer
 *
 * @deprecated Untested!!
 *
 * @todo Test!
 *
 * @todo Comment!
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
    
  void setType(int32_t t){type=t;}
  void setPosition(VectorInt3 pos){position=pos;}
  void setPosition(int32_t x,int32_t y,int32_t z){position.setX(x);position.setY(y);position.setZ(z);}
  int32_t getType() const{return type;}
  const VectorInt3 getPosition() const {return position;}
  size_t getParticleIndex() const {return particleIndex;}
  
private:
  VectorInt3 position;
  int32_t type;
  size_t particleIndex; /*set when apply is called. useful if Features want to alter the particle when applying the move*/
};



/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/*
 * @brief reset the probability
 */
template <class IngredientsType>
void MoveAddScMonomer::init(const IngredientsType& ing)
{
    this->resetProbability();
    particleIndex=ing.getMolecules().size();
}

/*****************************************************************************/
/*
 * check if the move is allowed by the system given as argument.
 */
template <class IngredientsType>
bool MoveAddScMonomer::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}
  
/*****************************************************************************/
/*
 * apply the move to the system given as argument
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

#endif
