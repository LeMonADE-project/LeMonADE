/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015,2021 by
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

#ifndef LEMONADE_UPDATER_MOVES_MOVECONNECTSCREACTIVE_H
#define LEMONADE_UPDATER_MOVES_MOVECONNECTSCREACTIVE_H
#include <limits>
#include <LeMonADE/updater/moves/MoveConnectBase.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveConnectScReactive
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 **/
/*****************************************************************************/

class MoveConnectScReactive:public MoveConnectBase<MoveConnectScReactive>
{
public:
  MoveConnectScReactive(){
      // old implementation until May2021
    /*shellPositions[0]=VectorInt3( 2, 0, 0);
    shellPositions[1]=VectorInt3(-2, 0, 0);
    shellPositions[2]=VectorInt3( 0, 2, 0);
    shellPositions[3]=VectorInt3( 0,-2, 0);
    shellPositions[4]=VectorInt3( 0, 0, 2);
    shellPositions[5]=VectorInt3( 0, 0,-2);
     */
    
    shellPositions[0]  = VectorInt3(2,0,0);
    shellPositions[1]  = VectorInt3(0,0,2);
    shellPositions[2]  = VectorInt3(0,2,0);
    shellPositions[3]  = VectorInt3(-2,0,0);
    shellPositions[4]  = VectorInt3(0,0,-2);
    shellPositions[5]  = VectorInt3(0,-2,0);

    shellPositions[6]  = VectorInt3(2,1,0);
    shellPositions[7]  = VectorInt3(1,0,2);
    shellPositions[8]  = VectorInt3(0,2,1);
    shellPositions[9]  = VectorInt3(-2,1,0);
    shellPositions[10] = VectorInt3(1,0,-2);
    shellPositions[11] = VectorInt3(0,-2,1);
    shellPositions[12] = VectorInt3(2,-1,0);
    shellPositions[13] = VectorInt3(-1,0,2);
    shellPositions[14] = VectorInt3(0,2,-1);
    shellPositions[15] = VectorInt3(-2,-1,0);
    shellPositions[16] = VectorInt3(-1,0,-2);
    shellPositions[17] = VectorInt3(0,-2,-1);
    shellPositions[18] = VectorInt3(1,2,0);
    shellPositions[19] = VectorInt3(2,0,1);
    shellPositions[20] = VectorInt3(0,1,2);
    shellPositions[21] = VectorInt3(-1,2,0);
    shellPositions[22] = VectorInt3(2,0,-1);
    shellPositions[23] = VectorInt3(0,-1,2);
    shellPositions[24] = VectorInt3(1,-2,0);
    shellPositions[25] = VectorInt3(-2,0,1);
    shellPositions[26] = VectorInt3(0,1,-2);
    shellPositions[27] = VectorInt3(-1,-2,0);
    shellPositions[28] = VectorInt3(-2,0,-1);
    shellPositions[29] = VectorInt3(0,-1,-2);

    shellPositions[30] = VectorInt3(2,1,1);
    shellPositions[31] = VectorInt3(1,1,2);
    shellPositions[32] = VectorInt3(1,2,1);
    shellPositions[33] = VectorInt3(-2,1,1);
    shellPositions[34] = VectorInt3(1,1,-2);
    shellPositions[35] = VectorInt3(1,-2,1);
    shellPositions[36] = VectorInt3(2,-1,1);
    shellPositions[37] = VectorInt3(-1,1,2);
    shellPositions[38] = VectorInt3(1,2,-1);
    shellPositions[39] = VectorInt3(-2,-1,1);
    shellPositions[40] = VectorInt3(-1,1,-2);
    shellPositions[41] = VectorInt3(1,-2,-1);
    shellPositions[42] = VectorInt3(2,1,-1);
    shellPositions[43] = VectorInt3(1,-1,2);
    shellPositions[44] = VectorInt3(-1,2,1);
    shellPositions[45] = VectorInt3(-2,1,-1);
    shellPositions[46] = VectorInt3(1,-1,-2);
    shellPositions[47] = VectorInt3(-1,-2,1);
    shellPositions[48] = VectorInt3(2,-1,-1);
    shellPositions[49] = VectorInt3(-1,-1,2);
    shellPositions[50] = VectorInt3(-1,2,-1);
    shellPositions[51] = VectorInt3(-2,-1,-1);
    shellPositions[52] = VectorInt3(-1,-1,-2);
    shellPositions[53] = VectorInt3(-1,-2,-1);
    
  }

  // overload initialise function to be able to set the moves index and direction if neccessary
  template <class IngredientsType> void init(const IngredientsType& ing);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, uint32_t partner);
  template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, VectorInt3 dir );

  template <class IngredientsType> bool check(IngredientsType& ing);
  template< class IngredientsType> void apply(IngredientsType& ing);

private:
  // holds the possible move directions
  /**
   * @brief Array that holds the 54 possible move directions
   *
   * @details 
   * * shellPositions   = (dx, dy, dz)
   * * shellPositions[..]=(+-2,  0,  0); //  6
   * * shellPositions[..]=(+-2,+-1,  0); // 24
   * * shellPositions[..]=(+-2,+-1,+-1); // 24
   */
  VectorInt3 shellPositions[54];
  const size_t shellPositionsEntries = 54;
  
  //! number of randomly choosen direction out of the shell
  const size_t numRandomSelectedDirections = 1;
  
  //! contains the id of the reactive partner in the randomly selected shell direction
  std::map<uint32_t, VectorInt3> IdReactivePartnerInShell;
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
void MoveConnectScReactive::init(const IngredientsType& ing)
{
  this->resetProbability();
  
  // clear the list of reactive partners in shell
  IdReactivePartnerInShell.clear();

  //draw index
  auto nUnreactedMonomers(ing.getNUnreactedMonomers());
  if( nUnreactedMonomers ==0 ) {
    this->setIndex(std::numeric_limits<uint32_t>::max());
    this->setPartner(std::numeric_limits<uint32_t>::max() );
  }
  else {
    // randomly select a monomer
    auto index(this->randomNumbers.r250_rand32() % nUnreactedMonomers);
    auto it (ing.getUnreactiveMonomers().cbegin());
    std::advance( it, index);
    this->setIndex(it->first);
    
    // search the vicinity for reactive monomers
    for(size_t test = 0; test < numRandomSelectedDirections; test++)
    {
        //draw direction
        VectorInt3 randomDir(shellPositions[ this->randomNumbers.r250_rand32() % shellPositionsEntries]);
        
        // Get the lattice value at a certain point - see FeatureConnectionSc
        // return monomer index between (0; molecules.size()-1) if there's a monomer
        // return uint32_t(-1)=4294967295 if place is empty
        uint32_t idPartner=ing.getIdFromLattice(ing.getMolecules()[this->getIndex()]+randomDir);
        
        // check if selected position contains a reactive partner with vaccant reaction site
        if ((std::numeric_limits<uint32_t>::max() != idPartner ) && (ing.getMolecules()[idPartner].isReactive()) && (ing.getMolecules().getNumLinks(idPartner) < ing.getMolecules()[idPartner].getNumMaxLinks()) )
        {
            IdReactivePartnerInShell[idPartner]=randomDir;
        }
    }
    
    // select a random partner and direction for connection move
    if(IdReactivePartnerInShell.size() == 0)
    {
        // no partner at all - move will rejected anyway
        this->setDir(shellPositions[0]); // it's arbitary 
        this->setPartner(std::numeric_limits<uint32_t>::max() );
    }
    else
    {   
        // select a random reactive monomer from the list
        auto idx(this->randomNumbers.r250_rand32() % IdReactivePartnerInShell.size());
        auto itID (IdReactivePartnerInShell.cbegin());
        std::advance( itID, idx);
        
        this->setDir(itID->second);
        this->setPartner(itID->first);
    }
    
  }
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
void MoveConnectScReactive::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();
  
  // clear the list of reactive partners in shell
  IdReactivePartnerInShell.clear();

  if ( !ing.checkCapableFormingBonds(index) )
  {
      this->setIndex( index );
      this->setDir(shellPositions[0]);
      this->setPartner(std::numeric_limits<uint32_t>::max());
  }else{
    //set index
    if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
      this->setIndex( index );
    else
      throw std::runtime_error("MoveConnectScReactive::init(ing, index): index out of range!");
      
    // search the vicinity for reactive monomers
    for(size_t test = 0; test < numRandomSelectedDirections; test++)
    {
        //draw direction
        VectorInt3 randomDir(shellPositions[ this->randomNumbers.r250_rand32() % shellPositionsEntries]);
        
        // Get the lattice value at a certain point - see FeatureConnectionSc
        // return monomer index between (0; molecules.size()-1) if there's a monomer
        // return uint32_t(-1)=4294967295 if place is empty
        uint32_t idPartner=ing.getIdFromLattice(ing.getMolecules()[this->getIndex()]+randomDir);
        
        // check if selected position contains a partner
        if ((std::numeric_limits<uint32_t>::max() != idPartner ) && ing.getMolecules()[idPartner].isReactive() )
        {
            IdReactivePartnerInShell[idPartner]=randomDir;
        }
    }
    
    // select a random partner and direction for connection move
    if(IdReactivePartnerInShell.size() == 0)
    {
        // no partner at all - move will rejected anyway
        this->setDir(shellPositions[0]); // it's arbitary 
        this->setPartner(std::numeric_limits<uint32_t>::max() );
    }
    else
    {   
        // select a random reactive monomer from the list
        auto idx(this->randomNumbers.r250_rand32() % IdReactivePartnerInShell.size());
        auto itID (IdReactivePartnerInShell.cbegin());
        std::advance( itID, idx);
        
        this->setDir(itID->second);
        this->setPartner(itID->first);
    }
  }
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
void MoveConnectScReactive::init(const IngredientsType& ing, uint32_t index, uint32_t partner)
{
  this->resetProbability();
  if ( !ing.checkCapableFormingBonds(index) )
  {
    this->setIndex( index );
    this->setDir(shellPositions[0]);
    this->setPartner(std::numeric_limits<uint32_t>::max());
  }else{
  
    //set index
    if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
      this->setIndex( index );
    else
      throw std::runtime_error("MoveConnectScReactive::init(ing, index): index out of range!");

    //draw direction
    this->setDir(ing.getMolecules()[partner].getVector3D()-ing.getMolecules()[index].getVector3D());
    this->setPartner(partner);
  }
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
void MoveConnectScReactive::init(const IngredientsType& ing, uint32_t index, VectorInt3 dir )
{
  this->resetProbability();
  if ( !ing.checkCapableFormingBonds(index) )
  {
      this->setIndex( index );
      this->setDir(shellPositions[0]);
      this->setPartner(std::numeric_limits<uint32_t>::max());
  }else{
    //set index
    if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
      this->setIndex( index );
    else
      throw std::runtime_error("MoveConnectScReactive::init(ing, index, bondpartner): index out of range!");

    //set direction if applicable
    if (std::find(std::begin(shellPositions), std::end(shellPositions), dir) != std::end(shellPositions))
    {
      this->setDir(dir);
    }
    else
      throw std::runtime_error("MoveLocalSc::init(ing, dir): direction vector out of range!");
    this->setPartner(ing.getIdFromLattice(ing.getMolecules()[this->getIndex()]+dir));
  }
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
bool MoveConnectScReactive::check(IngredientsType& ing)
{
  /**
   * @todo Think about this approach. Maybe we can put this statement somewhere else? 
   */
  if (std::numeric_limits<uint32_t>::max() == this->getPartner() ) return false ; 
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
void MoveConnectScReactive::apply(IngredientsType& ing)
{
	///@todo Think about the applying of move. Esp. make this independent of the order to avoid confusion!!
	///@todo check if it makes any difference in this case?!

	//FIRST bond is inserted 
	ing.modifyMolecules().connect(this->getIndex(),this->getPartner());
	
	//THEN all features are applied 
  	ing.applyMove(ing,*this);
	
	
	//the next line can produce a lot of output
	//thus use only if needed:
	//std::cout << " Connect: " << this->getIndex() << " to " << this->getPartner() << " at ("<<ing.modifyMolecules()[this->getIndex()] << ") - ("<< ing.modifyMolecules()[this->getPartner()] <<")"<< std::endl;
}

#endif
