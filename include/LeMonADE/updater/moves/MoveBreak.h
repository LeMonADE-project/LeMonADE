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

#ifndef LEMONADE_UPDATER_MOVES_MOVEBREAK_H
#define LEMONADE_UPDATER_MOVES_MOVEBREAK_H
#include <limits>
#include "./MoveBreakBase.h"

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveBreak
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/

class MoveBreak:public MoveBreakBase<MoveBreak>
{
public:
  MoveBreak(){};

  // overload initialise function to be able to set the moves index and direction if neccessary
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, uint32_t partner);

  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);

private:

};



/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template <class IngredientsType>
void MoveBreak::init(const IngredientsType& ing)
{
  this->resetProbability();

  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );
  //draw bond partner 
  auto bondID(this->randomNumbers.r250_rand32() %  ing.getMolecules().getNumLinks(this->getIndex()));
  this->setPartner( ing.getMolecules().getNeighborIdx(this->getIndex(), bondID) );
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 **/
template <class IngredientsType>
void MoveBreak::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index < ing.getMolecules().size() ) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveBreak::init(ing, index): index out of range!");

  //draw bond partner 
  auto bondID(this->randomNumbers.r250_rand32() %  ing.getMolecules().getNumLinks(this->getIndex()));
  this->setPartner( ing.getMolecules().getNeighborIdx(this->getIndex(), bondID) );
  
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 * @param partner index of the monomer to which index is connected
 **/
template <class IngredientsType>
void MoveBreak::init(const IngredientsType& ing, uint32_t index, uint32_t partner)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveBreak::init(ing, index, partner): index out of range!");

  //draw bond partner 
  if( (partner >= 0) && (partner <= (ing.getMolecules().size()-1)) )
  {
      this->setPartner(partner);
  }
  else
    throw std::runtime_error("MoveBreak::init(ing, index, partner): partner out of range!");
  
}
/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 **/
template <class IngredientsType>
bool MoveBreak::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}
/*****************************************************************************/
/**
 * @brief Apply the move to the system , e.g. add the displacement to Vertex (monomer) position.
 *
 * @details As first step: all Feature should apply the move using applyMove().\n
 * Second: Modify the positions etc. of the Vertex etc.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template< class IngredientsType>
void MoveBreak::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!
	///@todo check if it makes any difference in this case?! 
	//FIRST disconnect
	ing.modifyMolecules().disconnect(this->getIndex(),this->getPartner());	
	//THEN all features are applied 
  	ing.applyMove(ing,*this);
}

#endif
