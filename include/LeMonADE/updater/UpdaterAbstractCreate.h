/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2016 by
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
 * @details This abstract class provides the three basic functions to create systems: add a single monomer, add a connected monomer and move the system to find some free space.
 * This Updater requires FeatureAttributes.
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>

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

  //! function to add a standalone monomer
  bool addSingleMonomer(int32_t type=1);

  //! function to add a monomer to a parent monomer
  bool addMonomerToParent(uint32_t parent_id, int32_t type=1);

  //! function to add a monomer at a spezific place
  bool addMonomerAtPosition(VectorInt3 pos, int32_t type=1);

  //! function to add a new monomer inbetween two others
  bool addMonomerInsideConnectedPair(uint32_t indexA, uint32_t indexB, int32_t type=1);

  //! function to add a new monomer at two others
  bool addMonomerAtConnectedPair(uint32_t indexA, uint32_t indexB, int32_t type=1);
  
  //! function to add a ring around a chain monomer
  bool addRing(uint32_t parent_id, int32_t type=1, uint32_t NRingMonomers=5);
  
  //! function to move the whole system
  void moveSystem(int32_t nsteps);

  //! function to find groups of connected monomers
  void linearizeSystem();

  //! function to get a random bondvector of length 2
  VectorInt3 randomBondvector();

private:
  RandomNumberGenerators rng;

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
bool UpdaterAbstractCreate<IngredientsType>::addSingleMonomer(int32_t type){
  // set properties of add Monomer Move
  MoveAddMonomerSc<> addmove;
  addmove.init(ingredients);
  addmove.setTag(type);

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
bool UpdaterAbstractCreate<IngredientsType>::addMonomerToParent(uint32_t parent_id, int32_t type){
  // set properties of add Monomer Move
  MoveAddMonomerSc<> addmove;
  addmove.init(ingredients);
  addmove.setTag(type);

  int32_t counter(0);

  while(counter<10000){
    //try at most 30 random bondvectors to find a new monomer position
    for(uint i=0;i<30;i++){
      VectorInt3 bondvector(randomBondvector());
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
    moveSystem(2);
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
bool UpdaterAbstractCreate<IngredientsType>::addMonomerAtPosition(VectorInt3 position, int32_t type){
  MoveAddMonomerSc<> addmove;
  addmove.init(ingredients);
  addmove.setTag(type);

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
bool UpdaterAbstractCreate<IngredientsType>::addMonomerInsideConnectedPair(uint32_t indexA, uint32_t indexB, int32_t type){
  //first check if monomers are connected
  if( ! ingredients.getMolecules().areConnected(indexA,indexB))
    return false;

  MoveAddMonomerSc<> addmove;
  addmove.init(ingredients);
  addmove.setTag(type);

  int32_t counter(0);

  while(counter<10000){
    //try at most 30 random bondvectors to find a new monomer position
    for(uint i=0;i<30;i++){
      // get a random bondvector
      VectorInt3 bondvector(randomBondvector());
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
    moveSystem(2);
    counter++;
  }
  return false;
}
/******************************************************************************/
/**
 * @brief function to add a new monomer to two already existing ones such that all three monomers are connected.
 * @param indexA
 * @param indexB
 * @param type attribute tag of the new monomer
 * @return <b false> if position is not free, <b true> if move was applied
 */
template<class IngredientsType>
bool UpdaterAbstractCreate<IngredientsType>::addMonomerAtConnectedPair(uint32_t indexA, uint32_t indexB, int32_t type){
  //first check if monomers are connected
  if( ! ingredients.getMolecules().areConnected(indexA,indexB))
    return false;
  
  MoveAddMonomerSc<> addmove;
  addmove.init(ingredients);
  addmove.setTag(type);

  int32_t counter(0);

  while(counter<10000){
    //try at most 30 random bondvectors to find a new monomer position
    for(uint i=0;i<30;i++){
      // get a random bondvector
      VectorInt3 bondvector(randomBondvector());
      // set position of new monomer
      addmove.setPosition(ingredients.getMolecules()[indexA]+bondvector);

      // check new position (excluded volume, other features)
      if(addmove.check(ingredients)==true){
	//check the new bondvector between the new monomer and indexB
	VectorInt3 checkBV(addmove.getPosition()-ingredients.getMolecules()[indexB]);
	if( (checkBV.getLength() < 3) && (ingredients.getBondset().isValidStrongCheck(checkBV)) ){
	  addmove.apply(ingredients);
	  ingredients.modifyMolecules().connect( indexA, (ingredients.getMolecules().size()-1) );
	  ingredients.modifyMolecules().connect( indexB, (ingredients.getMolecules().size()-1) );
	  return true;
	}
      }
    }
    // if no position matches, we need to move the system a bit
    moveSystem(2);
    counter++;
  }
  return false;
}
/**
 * @brief add a ring threaded on a chain at position of parent monomer
 * 
 * @todo are there conformations where bonds of the ring cross the bonds of the chain?
 * 
 * @param parent 
 * @param type monomer tag 
 * @param NRingMonomers number of monomers per ring 
 */
template < class IngredientsType >
bool UpdaterAbstractCreate<IngredientsType>::addRing(uint32_t parent, int32_t type, uint32_t NRingMonomers)
  {
    if(parent < ingredients.getMolecules().size()) 
    {
        uint32_t attempts(0);
        while (attempts<10000)
        {
	  uint32_t neighborDirection(0);
	  uint32_t currentParent(parent);
	  while (true)
	  {
	    uint32_t neighborID(ingredients.getMolecules().getNeighborIdx(currentParent,neighborDirection));
	    if(ingredients.getMolecules().getNumLinks(neighborID) != 2 )
	    {
		if(neighborDirection==1) break;
		else {neighborDirection++; currentParent=parent;}
	    }
	    else
	    {
		VectorInt3 neighborPosition(ingredients.getMolecules()[neighborID]);
		VectorInt3 parentPosition(ingredients.getMolecules()[currentParent]);
		VectorInt3 bond(neighborPosition-parentPosition);
		uint32_t otherParentNeighbor; 
		uint32_t otherNeighborNeighbor;
		neighborDirection ? otherParentNeighbor=(ingredients.getMolecules().getNeighborIdx(currentParent,0)) :  otherParentNeighbor=(ingredients.getMolecules().getNeighborIdx(currentParent,1));
		(ingredients.getMolecules().getNeighborIdx(neighborID,0)==currentParent) ? otherNeighborNeighbor=(ingredients.getMolecules().getNeighborIdx(neighborID,1)): otherNeighborNeighbor=(ingredients.getMolecules().getNeighborIdx(neighborID,0));
		VectorInt3 neighborParentBond(parentPosition-ingredients.getMolecules()[otherParentNeighbor]);
		VectorInt3 neighborNeighborBond(neighborPosition-ingredients.getMolecules()[otherNeighborNeighbor]);
		if(bond.getLength() == 2 && (neighborParentBond.getLength()<3) &&(neighborNeighborBond.getLength()<3))
		{
		    for(uint32_t i=0; i<8; i++)
		    {
			
			VectorInt3 StartPosition=parentPosition;
			//set the possible conformations 
			int32_t dx1, dx2;
			int32_t dy1, dz1;
			int32_t dy2, dz2;
			VectorInt3 vec1,vec2;
			if(bond.getX() == 2 || bond.getX() == -2)
			{
			  StartPosition+=VectorInt3(bond.getX()/2,0,0);
			  dx1=0;dx2=0;//dy1=2;dy2=0;dz1=0;dz2=2; 
			  if(rng.r250_drand()>0.5){dy1=2;dy2=0;dz1=0;dz2=2;}
			  else {dy1=0;dy2=2;dz1=2;dz2=0;}
			  int32_t dx3,dx4;
			  if(rng.r250_drand()>0.5){dx3=1;dx4=-1;}
			  else{dx3=-1;dx4=1;}
// 			  rng.r250_drand()>0.5 ? dx3=1 : dx3=-1;
// 			  rng.r250_drand()>0.5 ? dx4=1 : dx4=-1;
			  vec1=VectorInt3( dx2+dx3,  dy2,  dz2);
			  vec2=VectorInt3(-dx2+dx4, -dy2, -dz2);
			}
			else if(bond.getY() == 2 || bond.getY() == -2)
			{
			  StartPosition+=VectorInt3(0,bond.getY()/2,0);
			  dy1=0;dy2=0; 
			  if(rng.r250_drand()>0.5){dx1=2;dx2=0;dz1=0;dz2=2;}
			  else {dx1=0;dx2=2;dz1=2;dz2=0;}
			  int32_t dy3,dy4;
			  if(rng.r250_drand()>0.5){dy3=1;dy4=-1;}
			  else{dy3=-1;dy4=1;}
// 			  rng.r250_drand()>0.5 ? dy3=1 : dy3=-1;
// 			  rng.r250_drand()>0.5 ? dy4=1 : dy4=-1;
			  vec1=VectorInt3( dx2,  dy2+dy3,  dz2);
			  vec2=VectorInt3(-dx2, -dy2+dy4, -dz2);
			}
			else if(bond.getZ() == 2 || bond.getZ() == -2)
			{
			  StartPosition+=VectorInt3(0,0,bond.getZ()/2);
			  dz1=0;dz2=0; 
			  if(rng.r250_drand()>0.5){dx1=2;dx2=0;dy1=0;dy2=2;}
			  else {dx1=0;dx2=2;dy1=2;dy2=0;}
			  int32_t dz3,dz4;
			  if(rng.r250_drand()>0.5){dz3=1;dz4=-1;}
			  else{dz3=-1;dz4=1;}
// 			  rng.r250_drand()>0.5 ? dz3=1 : dz3=-1;
// 			  rng.r250_drand()>0.5 ? dz4=1 : dz4=-1;
			  vec1=VectorInt3( dx2,  dy2,  dz2+dz3);
			  vec2=VectorInt3(-dx2, -dy2, -dz2+dz4);
			}
			// set positions of monomers and two which guarantee that ring is threaded 
			std::vector<VectorInt3> PotentialPositions(6,StartPosition);
			PotentialPositions[0]+=VectorInt3( dx1, dy1, dz1);
			PotentialPositions[2]+=VectorInt3(-dx1,-dy1,-dz1);
			PotentialPositions[4]+=VectorInt3( dx2, dy2, dz2);//position to check, not to use for a monomer 
			PotentialPositions[5]+=VectorInt3(-dx2,-dy2,-dz2);//position to check, not to use for a monomer 
			PotentialPositions[1]+=vec1;
			PotentialPositions[3]+=vec2;
			
			// check if positions are occupied
			bool PositionsFit(true);
			for(uint32_t i=0; i <PotentialPositions.size();i++)
			{
			  MoveAddMonomerSc<> addmove;
			  addmove.init(ingredients);
			  addmove.setPosition(PotentialPositions[i]);
			  if(addmove.check(ingredients)==false){PositionsFit=false;}
			}
			// add ring to system
			if (PositionsFit)
			{
			  bool AddRing(true);
			  for(uint32_t i=0; i<NRingMonomers;i++) 
			  {
			    if ( i<4  && AddRing ) if(!addMonomerAtPosition(PotentialPositions[i], type) ) {AddRing=false;throw std::runtime_error("Could not add 4-monomer ring!");}
			    if ( i==3 && AddRing ) for(uint32_t i=0; i < 4 ;i++) ingredients.modifyMolecules().connect(((ingredients.getMolecules().size()-1) -i%4),((ingredients.getMolecules().size()-1)-(1+i)%4));
			    if ( i>3  && AddRing ) 
			    {
			      for (uint32_t j=0;j<i;j++)
			      {
				//this is a dangerous thing, because it could happen that a monomer is added that freeze in 4 monomers of the system. 
				//I must get new safer function ....
				if (addMonomerInsideConnectedPair(ingredients.getMolecules().size()-(1+j%i),ingredients.getMolecules().size()-(1+(1+j)%i),type)) 
				{
				  AddRing=true;
				  break;
				}else { AddRing=false;}
			      }
			    }
			  }

			  if (AddRing) return true;
			  else throw std::runtime_error("Could not add a ring!");
			  
			}
		    }
		    
		    //if still here: try again at next position		
		    currentParent=neighborID;
			
		}
		else {currentParent=neighborID;}
	    }
	  }
	  attempts++;
	  moveSystem(2);
        }
        return false; 
    }else{
      std::stringstream error;
      error<<"Given parent ID "<<parent<<" does not exist!";
      throw std::runtime_error(error.str());}
  };

/******************************************************************************/
/**
 * @brief function to move the whole system
 * @param nsteps number of MCS to move
 */
template<class IngredientsType>
void UpdaterAbstractCreate<IngredientsType>::moveSystem(int32_t nsteps){
  MoveLocalSc move;
  for(int32_t n=0;n<nsteps;n++){
    for(int32_t m=0;m<ingredients.getMolecules().size();m++){
      move.init(ingredients);
      if(move.check(ingredients)==true){
	move.apply(ingredients);
      }
    }
  }
}

/******************************************************************************/
/**
 * @brief function to find groups of connected monomers and resort ingredients
 * to write out longest possible bondvector series
 */

template<class IngredientsType>
void UpdaterAbstractCreate<IngredientsType>::linearizeSystem(){
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

}

/******************************************************************************/
/**
 * @brief function to get a random bondvector of length 2 which is part of the bondvectorset
 * @param nsteps number of MCS to move
 */
template<class IngredientsType>
VectorInt3 UpdaterAbstractCreate<IngredientsType>::randomBondvector(){
  //get a random direction for the bondvector of length 2
  uint32_t randBondVectorID((rng.r250_rand32() % 6));
  VectorInt3 bondvector;
  switch (randBondVectorID)
  {
    case 0:
      bondvector.setAllCoordinates(2,0,0);
      break;
    case 1:
      bondvector.setAllCoordinates(0,2,0);
      break;
    case 2:
      bondvector.setAllCoordinates(0,0,2);
      break;
    case 3:
      bondvector.setAllCoordinates(-2,0,0);
      break;
    case 4:
      bondvector.setAllCoordinates(0,-2,0);
      break;
    case 5:
      bondvector.setAllCoordinates(0,0,-2);
      break;
  }

  //check if bondvector is part of the bondvectorset
  if(ingredients.getBondset().isValid(bondvector)){
    return bondvector;
  }else
    throw std::runtime_error("UpdaterAbstractCreate::randomBondvector: bondvectors not part of the bondvectorset.");
}
#endif /* LEMONADE_UPDATER_ABSTRACT_CREATE_H */
