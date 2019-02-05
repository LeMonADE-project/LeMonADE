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

#ifndef LEMONADE_UPDATER_MOVES_MOVELOCALSC_H
#define LEMONADE_UPDATER_MOVES_MOVELOCALSC_H

#include <LeMonADE/updater/moves/MoveConnectBase.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveConnectSc
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/
class MoveConnectSc:public MoveConnectBase<MoveConnectSc>
{
public:
  MoveConnectSc(){
    shellPositions[0]=VectorInt3( 2, 0, 0);
    shellPositions[1]=VectorInt3(-1, 0, 0);
    shellPositions[2]=VectorInt3( 0, 2, 0);
    shellPositions[3]=VectorInt3( 0,-1, 0);
    shellPositions[4]=VectorInt3( 0, 0, 2);
    shellPositions[5]=VectorInt3( 0, 0,-1);
	
// 	shellPositions[ 0]=VectorInt3(2,0,0);
// 	shellPositions[ 1]=VectorInt3(2,0,1);
// 	shellPositions[ 2]=VectorInt3(2,1,0);
// 	shellPositions[ 3]=VectorInt3(2,1,1);
// 	
// 	shellPositions[ 4]=VectorInt3(-1,0,0);
// 	shellPositions[ 5]=VectorInt3(-1,0,1);
// 	shellPositions[ 6]=VectorInt3(-1,1,0);
// 	shellPositions[ 7]=VectorInt3(-1,1,1);
// 	
// 	shellPositions[ 8]=VectorInt3(0,2,0);
// 	shellPositions[ 9]=VectorInt3(1,2,0);
// 	shellPositions[10]=VectorInt3(0,2,1);
// 	shellPositions[11]=VectorInt3(1,2,1);
// 	
// 	shellPositions[12]=VectorInt3(0,-1,0);
// 	shellPositions[13]=VectorInt3(1,-1,0);
// 	shellPositions[14]=VectorInt3(0,-1,1);
// 	shellPositions[15]=VectorInt3(1,-1,1);
// 	
// 	shellPositions[16]=VectorInt3(0,0,2);
// 	shellPositions[17]=VectorInt3(0,1,2);
// 	shellPositions[18]=VectorInt3(1,0,2);
// 	shellPositions[19]=VectorInt3(1,1,2);
// 	
// 	shellPositions[20]=VectorInt3(0,0,-1);
// 	shellPositions[21]=VectorInt3(0,1,-1);
// 	shellPositions[22]=VectorInt3(1,0,-1);
// 	shellPositions[23]=VectorInt3(1,1,-1);
  }

  // overload initialise function to be able to set the moves index and direction if neccessary
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, uint32_t bondpartner);

  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);

private:
  // holds the possible move directions
  /**
   * @brief Array that holds the 6 possible move directions
   *
   * @details In the scBFM where each monomer occupies 8 lattice sites and only one position is store (front, bottom, left) 
   * there are 1 positions in every move direction, which could be checked. But the positive directions have an offset due to the 
   * * shellPositions   = (dx, dy, dz)
   * * shellPositions[0]= ( 2,  0,  0);
   * * shellPositions[1]= (-1,  0,  0);
   * * shellPositions[2]= ( 0,  2,  0);
   * * shellPositions[3]= ( 0, -1,  0);
   * * shellPositions[4]= ( 0,  0,  2);
   * * shellPositions[5]= ( 0,  0, -1);
   * -1 {-1,0,1}
   * 0 {-1,0,1}
   * 1 {-1,0,1}
   */
  VectorInt3 shellPositions[6];
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
void MoveConnectSc::init(const IngredientsType& ing)
{
  this->resetProbability();

  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );

  //draw direction
  uint32_t randomDir=this->randomNumbers.r250_rand32() % 6;
  this->setPartner( ing.getLatticeEntry(ingredients.getMolecules()[Mon] + shellPositions[randomDir]) );
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
void MoveConnectSc::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveConnectSc::init(ing, index): index out of range!");

  //draw direction
  uint32_t randomDir=this->randomNumbers.r250_rand32() % 6;
  this->setPartner( ing.getLatticeEntry(ingredients.getMolecules()[Mon] + shellPositions[randomDir]) );
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 * @param bondpartner index of the monomer to connect to 
 **/
template <class IngredientsType>
void MoveConnectSc::init(const IngredientsType& ing, uint32_t index, uint32_t bondpartner )
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveConnectSc::init(ing, index, bondpartner): index out of range!");

  //set bond partner
  if( (bondpartner >= 0) && (bondpartner <= (ing.getMolecules().size()-1)) )
    this->setPartner( bondpartner );
  else
    throw std::runtime_error("MoveConnectSc::init(ing, index, bondpartner): bondpartner out of range!");
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
bool MoveConnectSc::check(IngredientsType& ing)
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
void MoveConnectSc::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!
	///@todo check if it makes any difference in this case?!

	//move must FIRST be applied to the features
	ing.applyMove(ing,*this);
	
	//THEN the bond is inserted 
	ing.modifyMolecules().connect(this->getIndex(),this->getPartner());
	

}

#endif
