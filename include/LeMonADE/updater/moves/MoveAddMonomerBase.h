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

#ifndef LEMONADE_UPDATER_MOVES_MOVEADDMOVEBASE_H
#define LEMONADE_UPDATER_MOVES_MOVEADDMOVEBASE_H

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveAddMonomerBase
 *
 * @brief Base class for all add monomer bfm-moves. Specialized versions are MoveAddScMonomer and MoveAddBccMonomer.
 *
 * @details This class provides a base add monomer moves. The implementation
 * details (i.e. MoveAddScMonomer, MoveAddBccMonomer, optimized versions,etc.) are given by the template parameter.
 * This way, using the "Curiously Recurring Template Pattern (CRTP) to avoid virtual functions
 * but still providing their functionality. A specialized local
 * move type derived from this class must then be constructed in the following
 * way: "class MySpecialMove:public MoveAddMonomerBase<MySpecialMove>".
 * It must implement the functions init,check and apply. Calls to the corresponding
 * functions in this base class are then redirected to the specialized implementation.
 * See for an example the class MoveAddScMonomer. 
 * Since this class simply serves as a common base, the functions apply(), check(), and init() don't do
 * anything particular.
 *
 * @tparam <SpecializedMove> name of the specialized move.
 *
 * @todo Test for this Base class?
 **/
/*****************************************************************************/
template <class SpecializedMove> 
class MoveAddMonomerBase:public MoveBase
{
public:
  //here come the functions that are implemented by the specialization
  //! Reset the probability
  template <class IngredientsType> void init(const IngredientsType& ingredients);
  //! Check if the move is allowed by the system given as argument.
  template <class IngredientsType> void check(const IngredientsType& ingredients);
  //! Apply the move to the system given as argument
  template <class IngredientsType> void apply(IngredientsType& ingredients);
  
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
  
protected:
  void setParticleIndex(size_t index){particleIndex=index;}
  
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

////////////////////////////////////////////////////////////////////////////////
// implementation of the members
////////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
/**
 * @brief Initialize the move. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType> 
void MoveAddMonomerBase<SpecializedMove>::init(const IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->init(ingredients);
}

/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter.
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType> 
void MoveAddMonomerBase<SpecializedMove>::check(const IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->check(ingredients);
}

/*****************************************************************************/
/**
 * @brief Apply the move to the system.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 * */
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType> 
void MoveAddMonomerBase<SpecializedMove>::apply(IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->apply(ingredients);
}

#endif