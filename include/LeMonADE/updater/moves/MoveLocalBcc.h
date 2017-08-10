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

#ifndef LEMONADE_UPDATER_MOVES_MOVELOCALBCC_H
#define LEMONADE_UPDATER_MOVES_MOVELOCALBCC_H

#include <LeMonADE/updater/moves/MoveLocalBase.h>


/*****************************************************************************/
/**
 * @file
 *
 * @class MoveLocalBcc
 *
 * @brief Standard local bfm-move on body-centered-cubic lattice for the bccBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 * In the bccBFM the 8 classic moves (dx,dy,dz) are along the diagonal of the lattice-axes as:
 * * movedir = (dx, dy, dz)
 * * movedir1 = ( 1,  1,  1)
 * * movedir2 = (-1,  1,  1)
 * * movedir3 = ( 1, -1,  1)
 * * movedir4 = ( 1,  1, -1)
 * * movedir5 = ( 1, -1, -1)
 * * movedir6 = (-1,  1, -1)
 * * movedir7 = (-1, -1,  1)
 * * movedir8 = (-1, -1, -1)
 *
 **/
/*****************************************************************************/

class MoveLocalBcc : public MoveLocalBase<MoveLocalBcc>
{
public:
	MoveLocalBcc(){};
	virtual ~MoveLocalBcc(){};

	template <class IngredientsType> void init(const IngredientsType& ing);
	template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
	template <class IngredientsType> void init(const IngredientsType& ing, VectorInt3 dir);
	template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, VectorInt3 dir);
	template <class IngredientsType> bool check(IngredientsType& ing);
	template< class IngredientsType> void apply(IngredientsType& ing);


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
 */
template <class IngredientsType>
void MoveLocalBcc::init(const IngredientsType& ing)
{
  this->resetProbability();

  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );

  //draw direction
  uint32_t randNumber=this->randomNumbers.r250_rand32();
  this->setDir(VectorInt3(
	((randNumber & 1) ? 1 :  -1),
	((randNumber & 2) ? 1 :  -1),
	((randNumber & 4) ? 1 :  -1)	)	);

}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * set the index of the moved monomer to the argument passed.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be moved
 */
template <class IngredientsType>
void MoveLocalBcc::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveLocalBcc::init(ing, index): index out of range!");

  //draw direction
  uint32_t randNumber=this->randomNumbers.r250_rand32();
  this->setDir(VectorInt3(
	((randNumber & 1) ? 1 :  -1),
	((randNumber & 2) ? 1 :  -1),
	((randNumber & 4) ? 1 :  -1)	)	);

}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given move direction.
 *
 * @details Resets the move probability to unity. Dice a Vertex (monomer) index
 * inside the graph and set the move direction to the argument passed.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param dir The direction of the move: must be one of the vectors P+-(1,1,1).
 */
template <class IngredientsType>
void MoveLocalBcc::init(const IngredientsType& ing, VectorInt3 dir)
{
  this->resetProbability();

  //draw index
  this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );

  //set direction
  if( (dir.getX()==1 || dir.getX()==-1 ) &&
    (dir.getY()==1 || dir.getY()==-1 ) &&
    (dir.getZ()==1 || dir.getZ()==-1 )  )
    this->setDir(dir);
  else
    throw std::runtime_error("MoveLocalBcc::init(ing, dir): direction vector out of range!");

}

/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 */
template <class IngredientsType>
void MoveLocalBcc::init(const IngredientsType& ing, uint32_t index, VectorInt3 dir)
{
  this->resetProbability();

  //draw index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveLocalBcc::init(ing, index, dir): index out of range!");

  //set direction
  if( (dir.getX()==1 || dir.getX()==-1 ) &&
    (dir.getY()==1 || dir.getY()==-1 ) &&
    (dir.getZ()==1 || dir.getZ()==-1 )  )
    this->setDir(dir);
  else
    throw std::runtime_error("MoveLocalBcc::init(ing, index, dir): direction vector out of range!");

}

/****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 *
 **/
template <class IngredientsType>
bool MoveLocalBcc::check(IngredientsType& ing)
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
void MoveLocalBcc::apply(IngredientsType& ing)
{
	///@todo make this independent of the order to avoid confusion!!

	//move must FIRST be applied to the features
	ing.applyMove(ing,*this);

	//THEN the position can be modified
	ing.modifyMolecules()[this->getIndex()]+=this->getDir();

}


#endif
