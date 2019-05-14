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
 * @brief Base class for all add monomer bfm-moves. Specialized versions are MoveAddMonomerSc and MoveAddMonomerBcc.
 *
 * @details This class provides a base add monomer moves. The implementation
 * details (i.e. MoveAddMonomerSc, MoveAddMonomerBcc, optimized versions,etc.) are given by the template parameter.
 * This way, using the "Curiously Recurring Template Pattern (CRTP) to avoid virtual functions
 * but still providing their functionality. A specialized local
 * move monomerTag derived from this class must then be constructed in the following
 * way: "class MySpecialMove:public MoveAddMonomerBase<MySpecialMove>".
 * It must implement the functions init,check and apply. Calls to the corresponding
 * functions in this base class are then redirected to the specialized implementation.
 * See for an example the class MoveAddMonomerSc.
 * Since this class simply serves as a common base, the functions apply(), check(), and init() don't do
 * anything particular.
 *
 * @tparam <SpecializedMove> name of the specialized move.
 *
 **/
/*****************************************************************************/
template<class SpecializedMove, class TagType>
class MoveAddMonomerBase:public MoveBase
{
public:
  MoveAddMonomerBase():monomerTag(TagType()),reactivity(false), numMaxLinks(0){};
  //here come the functions that are implemented by the specialization
  //! Reset the probability
  template <class IngredientsType> void init(const IngredientsType& ingredients);
  //! Check if the move is allowed by the system given as argument.
  template <class IngredientsType> void check(const IngredientsType& ingredients);
  //! Apply the move to the system given as argument
  template <class IngredientsType> void apply(IngredientsType& ingredients);

  //! setter function for monomerTag of new monomer
  void setTag(TagType tag){monomerTag=tag;}
  //! setter function for position of new monomer taking a VectorInt3
  void setPosition(VectorInt3 pos){position=pos;}
  //! setter function for position of new monomer taking a triple of ints
  void setPosition(int32_t x,int32_t y,int32_t z){position.setX(x);position.setY(y);position.setZ(z);}
    //! setter function for the reactivity 
  void setReactive(bool reactivity_){reactivity=reactivity_;}
  //! setter function for the maximum number of links the monomer should have 
  void setNumMaxLinks(uint32_t numMaxLinks_){numMaxLinks=numMaxLinks_;}
  //! getter function for the monomerTag of the new monomer
  TagType getTag() const{return monomerTag;}
  //! getter function for the position of the new monomer returning a VectorInt3
  const VectorInt3 getPosition() const {return position;}
  //! getter function for index of the new monomer. This is ing.getMolecules().size() before applying and ing.getMolecules().size()-1 after applying the move.
  size_t getMonomerIndex() const {return monomerIndex;}
  //! getter function for the reactivity 
  bool isReactive() const {return reactivity;}
  //! getter function for the maximum number of links the monomer should have 
  uint32_t getNumMaxLinks() const {return numMaxLinks;}
  
protected:
  void setMonomerIndex(size_t index){monomerIndex=index;}

private:
  //! position where the new monomer is placed in the simulation box
  VectorInt3 position;
  //! monomerTag that is applied to the new monomer and is required in  feature FeatureAttributes with TagType
  TagType monomerTag;
  /**
   * @brief Index of new PartileThis is ing.getMolecules().size() before applying and ing.getMolecules().size()-1 after applying the move.
   * @details It is set when apply is called. useful if Features want to alter the particle when applying the move
   */
  size_t monomerIndex;
  //! choose if the monomer is reactive ( can establish additional bonds in a connection process) used in feature FeatureConnectionSc
  bool reactivity;
  //! the maximum number of links the monomer should have mainly used in feature FeatureConnectionSc
  uint32_t numMaxLinks;
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
template<class SpecializedMove, class TagType>
template <class IngredientsType>
void MoveAddMonomerBase<SpecializedMove, TagType>::init(const IngredientsType& ingredients)
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
template<class SpecializedMove, class TagType>
template <class IngredientsType>
void MoveAddMonomerBase<SpecializedMove, TagType>::check(const IngredientsType& ingredients)
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
template<class SpecializedMove, class TagType>
template <class IngredientsType>
void MoveAddMonomerBase<SpecializedMove, TagType>::apply(IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->apply(ingredients);
}

#endif //LEMONADE_UPDATER_MOVES_MOVEADDMOVEBASE_H
