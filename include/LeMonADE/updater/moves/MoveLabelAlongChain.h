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

#ifndef LEMONADE_UPDATER_MOVES_MOVELABELALONGCHAING_H
#define LEMONADE_UPDATER_MOVES_MOVELABELALONGCHAING_H

#include "MoveLabelBase.h"

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveLabelAlongChain
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/
class MoveLabelAlongChain:public MoveLabelBase<MoveLabelAlongChain>
{
public:
  MoveLabelAlongChain(){}

  // overload initialise function to be able to set the moves index and direction if neccessary
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, int32_t dir);

  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);

private:
  // holds the possible move directions
  /**
   * @brief Array that holds the 6 possible move directions
   *
   * @details The label should stay along the chains and has the possibility to either go to the left(-1) or right(1).
   */
public:
  //! define an enum for the directions
  enum DIRECTIONS{
    LEFT=-1,       
    RIGHT=1,      
    UNDEFINED=0
  };
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
void MoveLabelAlongChain::init(const IngredientsType& ing)
{
  this->resetProbability();

  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );
  
  //draw direction
  int32_t randomDir=1-2*(this->randomNumbers.r250_rand32() % 2);
  this->setDir(randomDir);
  
  this->setConnectedLabel(ing.getLabelPartner(this->getIndex()));
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 **/
template <class IngredientsType>
void MoveLabelAlongChain::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveLabelAlongChain::init(ing, index): index out of range!");

  //draw direction
  int32_t randomDir=1-2*(this->randomNumbers.r250_rand32() % 2);
  this->setDir(randomDir);
  
  this->setConnectedLabel(ing.getLabelPartner(index));
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index and a given direction.
 *
 * @details Resets the move probability to unity and set the move properties index and direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 * @param dir Diretion given by and integer left(-1)-LEFT,right(1)-RIGHT.
 **/
template <class IngredientsType>
void MoveLabelAlongChain::init(const IngredientsType& ing, uint32_t index, int32_t dir)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveLabelAlongChain::init(ing, index, dir): index out of range!");

  //set direction
  if(dir==LEFT ||
     dir==RIGHT )
    this->setDir(dir);
  else
    throw std::runtime_error("MoveLabelAlongChain::init(ing, index, dir): direction is not valid!");
  
  this->setConnectedLabel(ing.getLabelPartner(index));
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
bool MoveLabelAlongChain::check(IngredientsType& ing)
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
void MoveLabelAlongChain::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!

	//move FIRST is applied to the features
	ing.applyMove(ing,*this);
	
	//if connected to another ring we need to check if the resultign new bond is valid
	if( this->getConnectedLabel() != 0 )
	{
	  ing.modifyMolecules().disconnect(this->getConnectedLabel()-1, this->getIndex()                );
	  ing.modifyMolecules().connect   (this->getConnectedLabel()-1, this->getIndex()+this->getDir() );
	}
	ing.modifyMolecules()[this->getIndex()+this->getDir()].setLabel(ing.getMolecules()[this->getIndex()].getLabel());
	ing.modifyMolecules()[this->getIndex()].setLabel(0);
}

#endif
