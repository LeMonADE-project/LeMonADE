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

#ifndef LEMONADE_UPDATER_MOVES_MOVELOCALSCDIAG_H
#define LEMONADE_UPDATER_MOVES_MOVELOCALSCDIAG_H

#include <LeMonADE/updater/moves/MoveLocalBase.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveLocalScDiag
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/
class MoveLocalScDiag:public MoveLocalBase<MoveLocalScDiag>
{
public:
  MoveLocalScDiag(){
	//perpendicular moves
	steps[0]=VectorInt3(1,0,0);
	steps[1]=VectorInt3(-1,0,0);
	steps[2]=VectorInt3(0,1,0);
	steps[3]=VectorInt3(0,-1,0);
	steps[4]=VectorInt3(0,0,1);
	steps[5]=VectorInt3(0,0,-1);
	//diagonal moves
	steps[6]=VectorInt3(0,1,1);
	steps[7]=VectorInt3(0,-1,1);
	steps[8]=VectorInt3(0,1,-1);
	steps[9]=VectorInt3(0,-1,-1);
	steps[10]=VectorInt3(1,0,1);
	steps[11]=VectorInt3(1,0,-1);
	steps[12]=VectorInt3(-1,0,1);
	steps[13]=VectorInt3(-1,0,-1);
	steps[14]=VectorInt3(1,1,0);
	steps[15]=VectorInt3(1,-1,0);
	steps[16]=VectorInt3(-1,1,0);
	steps[17]=VectorInt3(-1,-1,0);    
  }
  
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  template <class IngredientsType> void init(const IngredientsType& ing, VectorInt3 dir);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, VectorInt3 dir);
  
  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);
    
private:
  // holds the possible move directions
  /**
   * @brief Array that holds the 6 possible move directions
   *
   * @details In the scBFM the classic moves (dx,dy,dz) are along the lattice-axes as:
   * * steps   = (dx, dy, dz)
   * * steps[0]= ( 1,  0,  0);
   * * steps[1]= (-1,  0,  0);
   * * steps[2]= ( 0,  1,  0);
   * * steps[3]= ( 0, -1,  0);
   * * steps[4]= ( 0,  0,  1);
   * * steps[5]= ( 0,  0, -1);
   * * Plus the moves diagonal with square length 2.  
   */
  VectorInt3 steps[18];
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
void MoveLocalScDiag::init(const IngredientsType& ing)
{
  this->resetProbability();
  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );
  //draw direction
  uint32_t randomDir=this->randomNumbers.r250_rand32() % 18;
  
  this->setDir(steps[randomDir]);

}

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 **/
template <class IngredientsType>
void MoveLocalScDiag::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();
  
  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveLocalScDiag::init(ing, index): index out of range!");
  
  //draw direction
  uint32_t randomDir=this->randomNumbers.r250_rand32() % 18;
  
  this->setDir(steps[randomDir]);

}

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param dir The direction of the move: must be one of the vectors P+-(1,0,0).
 **/
template <class IngredientsType>
void MoveLocalScDiag::init(const IngredientsType& ing, VectorInt3 dir)
{
  this->resetProbability();
  
  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );
  
 //set direction
  if(dir==steps[ 0] ||
     dir==steps[ 1] ||
     dir==steps[ 2] ||
     dir==steps[ 3] ||
     dir==steps[ 4] ||
     dir==steps[ 5] ||
     dir==steps[ 6] ||
     dir==steps[ 7] ||
     dir==steps[ 8] ||
     dir==steps[ 9] ||
     dir==steps[10] ||
     dir==steps[11] ||
     dir==steps[12] ||
     dir==steps[13] ||
     dir==steps[14] ||
     dir==steps[15] ||
     dir==steps[16] ||
     dir==steps[17] 
    )
    this->setDir(dir);
  else
    throw std::runtime_error("MoveLocalScDiag::init(ing, dir): direction vector out of range!");
}

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 * @param dir The direction of the move: must be one of the vectors P+-(1,0,0).
 **/
template <class IngredientsType>
void MoveLocalScDiag::init(const IngredientsType& ing, uint32_t index, VectorInt3 dir)
{
  this->resetProbability();
  
  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveLocalScDiag::init(ing, index, dir): index out of range!");
  
 //set direction
  if(dir==steps[ 0] ||
     dir==steps[ 1] ||
     dir==steps[ 2] ||
     dir==steps[ 3] ||
     dir==steps[ 4] ||
     dir==steps[ 5] ||
     dir==steps[ 6] ||
     dir==steps[ 7] ||
     dir==steps[ 8] ||
     dir==steps[ 9] ||
     dir==steps[10] ||
     dir==steps[11] ||
     dir==steps[12] ||
     dir==steps[13] ||
     dir==steps[14] ||
     dir==steps[15] ||
     dir==steps[16] ||
     dir==steps[17] 
    )
    this->setDir(dir);
  else
    throw std::runtime_error("MoveLocalScDiag::init(ing, index, dir): direction vector out of range!");
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
bool MoveLocalScDiag::check(IngredientsType& ing)
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
void MoveLocalScDiag::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!
	
	//move must FIRST be applied to the features
	ing.applyMove(ing,*this);	

	//THEN the position can be modified
	ing.modifyMolecules()[this->getIndex()]+=this->getDir();
  
}

#endif
