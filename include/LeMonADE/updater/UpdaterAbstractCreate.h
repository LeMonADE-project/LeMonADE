/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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
#ifndef LEMONADE_UPDATER_ABSTRACT_CREATE_H
#define LEMONADE_UPDATER_ABSTRACT_CREATE_H


/**
 * @file
 *
 * @class UpdaterAbstractCreate
 *
 * @brief abstract updater to create systems
 * 
 * @details This abstract class provides the three basic functions to create systems: add a single monomer, add a connected monomer and move the system to find some free space
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/Vector3D.h>

template<class IngredientsType>
class UpdaterAbstractCreate:public AbstractUpdater
{
public:
  UpdaterAbstractCreate(IngredientsType& ingredients_):ingredients(ingredients_){}
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
protected:
  IngredientsType& ingredients;
  RandomNumberGenerators rng;
  
  //! function to add a monomer to a parent monomer
  bool add_monomer_to_parent(uint32_t parent_id, int32_t type);
  //! function to add a standalone monomer
  bool add_lonely_monomer(int32_t type);
  //! function to move the whole system
  void move_system(int32_t nsteps);
  //! function to add a monomer at a spezific place
  bool add_monomer_to_position(VectorInt3 pos, int32_t type);
  
};

/**
* The initialize function handles the new systems information.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAbstractCreate<IngredientsType>::initialize(){
  
}

/**
* Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterAbstractCreate<IngredientsType>::execute(){
  
}

/**
* Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAbstractCreate<IngredientsType>::cleanup(){
  
}

/******************************************************************************/
/**
 * @brief function to add a monomer to a parent monomer
 * @param parent_id id of monomer to connect with
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_monomer_to_parent(uint32_t parent_id, int32_t type=1){
  // set properties of add Monomer Move
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  
  int32_t counter(0);
  
  while(counter<10000){
    //try at most 30 random bondvectors to find a new monomer position
    for(uint i=0;i<30;i++){
      // get id of random bondvector with index <= 22, also P(2,0,0)
      std::size_t randBondVectorID((rng.r250_rand32() % 6)+17);
      VectorInt3 bondvector(ingredients.getBondset().getBondVector(randBondVectorID));
      // set position of new monomer
      addmove.setPosition(ingredients.getMolecules()[parent_id]+bondvector);
    
      // check new position (excluded volume)
      if(addmove.check(ingredients)==true){
	addmove.setType(type);
	addmove.apply(ingredients);
	ingredients.modifyMolecules().connect( parent_id, (ingredients.getMolecules().size()-1) );
	return true;
      }
    }
    // if no position matches, we need to move the system a bit
    move_system(2);
    counter++;
  }
  return false;
}

/******************************************************************************/
/**
 * @brief function to add a standalone monomer
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_lonely_monomer(int32_t type=1){
  // set properties of add Monomer Move
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  
  int32_t counter(0);
  while(counter<10000){
    VectorInt3 newPosition((rng.r250_rand32() % (ingredients.getBoxX()-1)),
				  (rng.r250_rand32() % (ingredients.getBoxY()-1)),
				  (rng.r250_rand32() % (ingredients.getBoxZ()-1)));
    addmove.setPosition(newPosition);
    if(addmove.check(ingredients)==true){
      addmove.setType(type);
      addmove.apply(ingredients);
      return true;
    }
    counter++;
  }
  return false;
}


/******************************************************************************/
/**
 * @brief function to move the whole system
 * @param nsteps number of MCS to move
 */
template<class IngredientsType>
void UpdaterAbstractCreate<IngredientsType>::move_system(int32_t nsteps){
  MoveLocalSc move;
  for(int32_t n=0;n<nsteps;n++){
    for(int32_t m=0;m<ingredients.getMolecules().size();m++){
      move.init(ingredients);
      if(move.check(ingredients)==true){
	move.apply(ingredients);
      }
    }
  }
  ingredients.modifyMolecules().setAge(ingredients.getMolecules().getAge()+nsteps);
}

/******************************************************************************/
/**
 * @brief function to add a monomer to a specific position. if position is not free return false 
 * @param VectorInt3 position
 * @return <b false> if position is not free, <b true> if move was applied
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_monomer_to_position(VectorInt3 position, int32_t type=1){
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  addmove.setPosition(position);
  if(addmove.check(ingredients)==true){
    addmove.setType(type);
    addmove.apply(ingredients);
    return true;
  }else{
    return false;
  }
  
}

/******************************************************************************/
/**
 * @brief function to rearrange monomers to get linear strands
 */
/*
template<class IngredientsType>
void UpdaterAbstractCreate<IngredientsType>::linearize_system(){
  //call ingredients copy constructor
  IngredientsType newIngredients(ingredients);
  
  //delete all the informations in molecules
  newIngredients.modifyMolecules().clear();
  
  //create a death iterator object, only one template parameter gives the default "always true"
  GraphIteratorDepthFirst<MoleculesType> iterator(ingredients.getMolecules());
  
  typedef std::vector < MonomerGroup<typename IngredientsType::molecules_type> > MonomerGroupVector;

  MonomerGroupVector LinearMonomerGroupsVector;

  fill_connected_groups( ingredients.getMolecules(), LinearMonomerGroupsVector, MonomerGroup<typename IngredientsType::molecules_type>(&(ingredients.getMolecules())), alwaysTrue );
  
  
}/* */

#endif /* LEMONADE_UPDATER_ABSTRACT_CREATE_H */
