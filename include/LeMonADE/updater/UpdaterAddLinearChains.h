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

#ifndef LEMONADE_UPDATER_SETUP_LINEAR_CHAINS
#define LEMONADE_UPDATER_SETUP_LINEAR_CHAINS
/**
 * @file
 *
 * @class UpdaterAddLinearChains
 *
 * @brief Updater to create a solution of monomdisperse linear chains.
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients
 * or a system with some monomers inside. This updater requires FeatureAttributes.
 * Two tags are added to the monomers in alternating manner, usually needed for GPU computing.
 *
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding eigther an empty simulation box for system setup
 * or a prefilled ingredients where the linear chains shall be added
 * @param NChain_ number of chains that are added to ingredients
 * @param NMonoPerChain_ number of monomer is each chain
 * @param type1_ attribute tag of "even" monomers
 * @param type2_ attribute tag of "odd" monomers
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>
#include <cmath>

template<class IngredientsType>
class UpdaterAddLinearChains: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
  UpdaterAddLinearChains(IngredientsType& ingredients_, uint32_t NChain_, uint32_t NMonoPerChain_, int32_t type1_=1, int32_t type2_=2, bool IsSolvent=false);

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

  //! getter function for write compressed solvent bool
  const bool getIsSolvent() const {return IsSolvent;}

  //! getter function for number of monomers in chains
  const int32_t getNMonomerPerChain() const {return NMonoPerChain;}

  //! getter function for number of chains
  const int32_t getNChain() const {return NChain;}

  //! getter function for calculated density
  const double getDensity() const {return density;}

private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;

  //! number of monomers in a chain
  uint32_t NMonoPerChain;

  //! number of linear chains in the box
  uint32_t NChain;

  //! lattice occupation density
  double density;

  //! bool for execution
  bool wasExecuted;

  //! attribute tag of even monomers
  int32_t type1;

  //! getAttributeTag of odd monomers
  int32_t type2;

  //! bool to check if chains of size 1 should be compressed to solvent
  bool IsSolvent;
};

/**
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param NChain_ number of chains to be added in the system instead of solvent
* @param NMonoPerChain_ number of monomers in one chain
*/
template < class IngredientsType >
UpdaterAddLinearChains<IngredientsType>::UpdaterAddLinearChains(IngredientsType& ingredients_, uint32_t NChain_, uint32_t NMonoPerChain_, int32_t type1_, int32_t type2_, bool IsSolvent_):
BaseClass(ingredients_), NChain(NChain_), NMonoPerChain(NMonoPerChain_), density(0.0), wasExecuted(false),
type1(type1_), type2(type2_), IsSolvent(IsSolvent_)
{}

/**
* @brief initialise function, calculate the target density to compare with at the end.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddLinearChains<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterAddLinearChains" << std::endl;

  // get the target density from the sum of existing monomers and the new added chains
  density=(double)( ingredients.getMolecules().size() + NMonoPerChain*NChain ) * 8  /(double)( ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ() );

  std::cout << "add "<<NChain*NMonoPerChain<<" monomers to the box"<<std::endl;

  execute();
}

/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterAddLinearChains<IngredientsType>::execute(){
  if(wasExecuted)
    return true;

  std::cout << "execute UpdaterAddLinearChains" << std::endl;

  //loop over chains and chain monomers and build it up
  for(uint32_t i=0;i<(NChain);i++){
    for(uint32_t j=0;j<(NMonoPerChain);j++){
      if(j==0)
	addSingleMonomer(type1);
      else{
	if(ingredients.getMolecules()[ingredients.getMolecules().size()-1].getAttributeTag() == type1)
	  addMonomerToParent(ingredients.getMolecules().size()-1,type2);
	else
	  addMonomerToParent(ingredients.getMolecules().size()-1,type1);
      }
    }
  }

  ingredients.synchronize();
  double lattice_volume(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());
  if(std::abs(density - ( (double)(ingredients.getMolecules().size()*8) / lattice_volume )) > 0.0000000001 ){
    std::cout << density << " " <<( (ingredients.getMolecules().size()*8) / lattice_volume)<<std::endl;
    throw std::runtime_error("UpdaterAddLinearChains: number of monomers in molecules does not match the calculated number of monomers!");
  }else{
    std::cout << "real lattice occupation density =" << (8*ingredients.getMolecules().size()) / lattice_volume<<std::endl;
    wasExecuted=true;
    // if we allow for colvent compression AND added single monomers (solvent): compress solvent
    if(IsSolvent && NMonoPerChain==1){
      ingredients.setCompressedOutputIndices(ingredients.getMolecules().size()-(NMonoPerChain*NChain),ingredients.getMolecules().size()-1);
    }else{
      linearizeSystem();
    }
    return true;
  }
}

/**
* @brief Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddLinearChains<IngredientsType>::cleanup(){

}


#endif /* LEMONADE_UPDATER_SETUP_LINEAR_CHAINS */