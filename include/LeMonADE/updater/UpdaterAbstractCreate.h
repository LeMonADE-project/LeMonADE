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
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DepthIterator.h>
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
  
  //! function to add a standalone monomer
  bool add_lonely_monomer(int32_t type=1);
  
  //! function to add a monomer to a parent monomer
  bool add_monomer_to_parent(uint32_t parent_id, int32_t type=1);
  
  //! function to add a monomer at a spezific place
  bool add_monomer_to_position(VectorInt3 pos, int32_t type=1);
  
  //! function to add a new monomer inbetween two others
  bool add_monomer_inside_connected_pair(uint32_t indexA, uint32_t indexB, int32_t type=1);
  
  //! function to move the whole system
  void move_system(int32_t nsteps);
  
  //! function to find groups of connected monomers
  void linearize_system();
  
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
 * @brief function to add a standalone monomer
 * @param type attribute tag of the new monomer
 * @return <b false> if position is not free, <b true> if move was applied
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_lonely_monomer(int32_t type){
  // set properties of add Monomer Move
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  addmove.setType(type);
  
  int32_t counter(0);
  while(counter<10000){
    VectorInt3 newPosition((rng.r250_rand32() % (ingredients.getBoxX()-1)),
				  (rng.r250_rand32() % (ingredients.getBoxY()-1)),
				  (rng.r250_rand32() % (ingredients.getBoxZ()-1)));
    addmove.setPosition(newPosition);
    if(addmove.check(ingredients)==true){
      addmove.apply(ingredients);
      return true;
    }
    counter++;
  }
  return false;
}

/******************************************************************************/
/**
 * @brief function to add a monomer to a parent monomer
 * @param parent_id id of monomer to connect with
 * @param type attribute tag of the new monomer
 * @return <b false> if position is not free, <b true> if move was applied
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_monomer_to_parent(uint32_t parent_id, int32_t type){
  // set properties of add Monomer Move
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  addmove.setType(type);
  
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
 * @brief function to add a monomer to a specific position. if position is not free return false 
 * @param VectorInt3 position
 * @param type attribute tag of the new monomer
 * @return <b false> if position is not free, <b true> if move was applied
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_monomer_to_position(VectorInt3 position, int32_t type){
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  addmove.setType(type);
  
  addmove.setPosition(position);
  if(addmove.check(ingredients)==true){
    addmove.apply(ingredients);
    return true;
  }else{
    return false;
  }
  
}

/******************************************************************************/
/**
 * @brief function to add a new monomer between two already existing ones instead of a bond.
 * @param indexA
 * @param indexB
 * @param type attribute tag of the new monomer
 * @return <b false> if position is not free, <b true> if move was applied
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::add_monomer_inside_connected_pair(uint32_t indexA, uint32_t indexB, int32_t type){
  //first check if monomers are connected
  if( ! ingredients.getMolecules().areConnected(indexA,indexB))
    return false;
  
  MoveAddScMonomer addmove;
  addmove.init(ingredients);
  addmove.setType(type);
  
  int32_t counter(0);
  
  while(counter<10000){
    //try at most 30 random bondvectors to find a new monomer position
    for(uint i=0;i<30;i++){
      // get id of random bondvector with index <= 22, also P(2,0,0)
      std::size_t randBondVectorID((rng.r250_rand32() % 6)+17);
      VectorInt3 bondvector(ingredients.getBondset().getBondVector(randBondVectorID));
      // set position of new monomer
      addmove.setPosition(ingredients.getMolecules()[indexA]+bondvector);
    
      // check new position (excluded volume, other features)
      if(addmove.check(ingredients)==true){
	//check the new bondvector bewten the new monomer and indexB
	VectorInt3 checkBV(addmove.getPosition()-ingredients.getMolecules()[indexB]);
	if( (checkBV.getLength() < 3) && (ingredients.getBondset().isValidStrongCheck(checkBV)) ){
	  addmove.apply(ingredients);
	  ingredients.modifyMolecules().connect( indexA, (ingredients.getMolecules().size()-1) );
	  ingredients.modifyMolecules().connect( indexB, (ingredients.getMolecules().size()-1) );
	  ingredients.modifyMolecules().disconnect( indexA, indexB );
	  return true;
	}
      }
    }
    // if no position matches, we need to move the system a bit
    move_system(2);
    counter++;
  }
  // create a meaningful return message
  std::cout << "UpdaterAbstractCreate::add_monomer_inside_connected_pair:  could not add new monomer between "<<indexA<< " and "<<indexB<<std::endl; 
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
 * @brief function to find groups of connected monomers and resort ingredients 
 * to write out longest possible bondvector series
 */

template<class IngredientsType>
void UpdaterAbstractCreate<IngredientsType>::linearize_system(){
  //call ingredients copy constructor
  IngredientsType newIngredients(ingredients);
  
  //delete all the informations in molecules
  newIngredients.modifyMolecules().clear();
  
  std::vector < MonomerGroup<typename IngredientsType::molecules_type> > LinearMonomerGroupsVector;

  fill_connected_groups( ingredients.getMolecules(), LinearMonomerGroupsVector, MonomerGroup<typename IngredientsType::molecules_type>(ingredients.getMolecules()), alwaysTrue() );
  
  for(size_t groups=0; groups < LinearMonomerGroupsVector.size(); ++groups){
    if(groups==0)
      newIngredients.modifyMolecules() = LinearMonomerGroupsVector[groups].copyGroup();
    else
      newIngredients.modifyMolecules() += LinearMonomerGroupsVector[groups].copyGroup();
  }
  
  ingredients=newIngredients;
  
}/* */

#endif /* LEMONADE_UPDATER_ABSTRACT_CREATE_H */
