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

#include "gtest/gtest.h"

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureExcludedVolumeScIdOnLattice.h>
#include <LeMonADE/utility/LatticePredicates.h>


using namespace std;

class TestFeatureExcludedVolumeScIdOnLattice: public ::testing::Test{
public:
  //dummy move class used to check response to unknown move type
  class UnknownMove:public MoveBase
  {
	  public:
		template<class IngredientsType> bool check(IngredientsType& ingredients) const
		{
		return ingredients.checkMove(ingredients,*this);
		}

		template<class IngredientsType> void apply(IngredientsType& ingredients)
		{
		ingredients.applyMove(ingredients,*this);
		}

		template <class IngredientsType> void init(const IngredientsType& ingredients){};
  };

  //redirect cout output
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };

private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
};

TEST_F(TestFeatureExcludedVolumeScIdOnLattice,LatticeEntryMonomerID)
{
  typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, MonomerID > MyExcludedVolumeFeature;
  typedef LOKI_TYPELIST_2(MyExcludedVolumeFeature, FeatureBondset<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;
  Ing ingredients;

  //prepare ingredients
    ingredients.setBoxX(32);
    ingredients.setBoxY(32);
    ingredients.setBoxZ(32);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);

    //one move of every type
    MoveLocalSc scmove;
    MoveAddMonomerSc<> addmove;
    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(1,0,30);
    
    // **************   synchronize lattice with monomer ID   **************
    ingredients.synchronize(ingredients);
    EXPECT_NO_THROW(ingredients.synchronize(ingredients));
    std::vector<VectorInt3> Position;
    Position.push_back(VectorInt3(0,0,0));
    Position.push_back(VectorInt3(2,0,0));
    Position.push_back(VectorInt3(1,0,30));
    
    //first monomer 
    VectorInt3 vec(Position[0] + VectorInt3(0,0,0));
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,1,0);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,0,1);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,-1,0);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
    //second monomer 
    vec=(Position[1] + VectorInt3(0,0,0));
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,1,0);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,0,1);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,-1,0);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(2, ingredients.getLatticeEntry(vec));
    //third monomer 
    vec=(Position[2] + VectorInt3(0,0,0));
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,1,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,0,1);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,-1,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    
    // **************   check sc move   **************
    scmove.init(ingredients,2,VectorInt3(0,1,0));
    scmove.check(ingredients);
    scmove.apply(ingredients);
    
    //third monomer: check if indecies are correctly moved
    vec=(Position[2]+ VectorInt3(0,1,0) + VectorInt3(0,0,0));
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,1,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,0,1);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,-1,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(3, ingredients.getLatticeEntry(vec));
    
    // **************   check sc add move   **************
    
    addmove.init(ingredients);
    
    addmove.setPosition(1,1,30);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(1,2,30);
    EXPECT_FALSE(addmove.check(ingredients));
    
    addmove.setPosition(1,3,30);
    EXPECT_TRUE(addmove.check(ingredients));
    
    //apply move and check position again
    addmove.apply(ingredients);
    addmove.init(ingredients);
    EXPECT_FALSE(addmove.check(ingredients));
    
    //fourth monomer: check if indecies are correctly moved
    vec=(VectorInt3(1,3,30) + VectorInt3(0,0,0));
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,1,0);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,0,1);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(1,0,0);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(0,-1,0);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    vec+=VectorInt3(-1,0,0);
    EXPECT_EQ(4, ingredients.getLatticeEntry(vec));
    
}

TEST_F(TestFeatureExcludedVolumeScIdOnLattice,LatticeEntryAttributeANDMonomerID)
{
	typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, AttributeANDMonomerID > MyExcludedVolumeFeature;
	typedef LOKI_TYPELIST_1(MyExcludedVolumeFeature) Features;

	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing ingredients;

	 //prepare ingredients1
	    ingredients.setBoxX(32);
	    ingredients.setBoxY(32);
	    ingredients.setBoxZ(32);
	    ingredients.setPeriodicX(1);
	    ingredients.setPeriodicY(1);
	    ingredients.setPeriodicZ(1);
	    //set up ingredients
	    typename Ing::molecules_type& molecules=ingredients.modifyMolecules();

	    molecules.resize(4);
	    molecules[0].setAllCoordinates(9,10,10);
	    molecules[0].setAttributeTag(0);
	    molecules[1].setAllCoordinates(7,8,12);
	    molecules[1].setAttributeTag(17);
	    molecules[2].setAllCoordinates(2,3,3);
	    molecules[2].setAttributeTag(28);
	    molecules[3].setAllCoordinates(30,17,13);
	    molecules[3].setAttributeTag(63);
	    ingredients.synchronize(ingredients);
	    EXPECT_NO_THROW(ingredients.synchronize(ingredients));
	    VectorInt3 vec(9,10,10);
	    EXPECT_EQ(1, ingredients.getLatticeEntry(vec)>>6);
	    EXPECT_EQ(0, ingredients.getLatticeEntry(vec) & 63);
	    vec.setAllCoordinates(7,8,12);
	    EXPECT_EQ(2, ingredients.getLatticeEntry(vec)>>6);
	    EXPECT_EQ(17, ingredients.getLatticeEntry(vec) & 63);
	    vec.setAllCoordinates(2,3,3);
	    EXPECT_EQ(3, ingredients.getLatticeEntry(vec)>>6);
	    EXPECT_EQ(28, ingredients.getLatticeEntry(vec) & 63);
	    vec.setAllCoordinates(30,17,13);
	    EXPECT_EQ(4, ingredients.getLatticeEntry(vec)>>6);
	    EXPECT_EQ(63, ingredients.getLatticeEntry(vec) & 63);
}


TEST_F(TestFeatureExcludedVolumeScIdOnLattice,LatticeEntryAttribute)
{
	typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, Attribute > MyExcludedVolumeFeature;
	typedef LOKI_TYPELIST_1(MyExcludedVolumeFeature) Features;

	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing ingredients;

	 //prepare ingredients1
	    ingredients.setBoxX(32);
	    ingredients.setBoxY(32);
	    ingredients.setBoxZ(32);
	    ingredients.setPeriodicX(1);
	    ingredients.setPeriodicY(1);
	    ingredients.setPeriodicZ(1);
	    //set up ingredients
	    typename Ing::molecules_type& molecules=ingredients.modifyMolecules();

	    molecules.resize(4);
	    molecules[0].setAllCoordinates(9,10,10);
	    molecules[0].setAttributeTag(0);
	    molecules[1].setAllCoordinates(7,8,12);
	    molecules[1].setAttributeTag(17);
	    molecules[2].setAllCoordinates(2,3,3);
	    molecules[2].setAttributeTag(28);
	    molecules[3].setAllCoordinates(30,17,13);
	    molecules[3].setAttributeTag(63);
	    ingredients.synchronize(ingredients);
	    EXPECT_NO_THROW(ingredients.synchronize(ingredients));
	    VectorInt3 vec(9,10,10);
	    EXPECT_EQ(0, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(7,8,12);
	    EXPECT_EQ(17, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(2,3,3);
	    EXPECT_EQ(28, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(30,17,13);
	    EXPECT_EQ(63, ingredients.getLatticeEntry(vec));
}

TEST_F(TestFeatureExcludedVolumeScIdOnLattice,Bool)
{
	typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, Bool > MyExcludedVolumeFeature;
	typedef LOKI_TYPELIST_1(MyExcludedVolumeFeature) Features;

	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing ingredients;

	 //prepare ingredients1
	    ingredients.setBoxX(32);
	    ingredients.setBoxY(32);
	    ingredients.setBoxZ(32);
	    ingredients.setPeriodicX(1);
	    ingredients.setPeriodicY(1);
	    ingredients.setPeriodicZ(1);
	    //set up ingredients
	    typename Ing::molecules_type& molecules=ingredients.modifyMolecules();

	    molecules.resize(4);
	    molecules[0].setAllCoordinates(9,10,10);
	    molecules[0].setAttributeTag(0);
	    molecules[1].setAllCoordinates(7,8,12);
	    molecules[1].setAttributeTag(17);
	    molecules[2].setAllCoordinates(2,3,3);
	    molecules[2].setAttributeTag(28);
	    molecules[3].setAllCoordinates(30,17,13);
	    molecules[3].setAttributeTag(63);
	    ingredients.synchronize(ingredients);
	    EXPECT_NO_THROW(ingredients.synchronize(ingredients));
	    VectorInt3 vec(9,10,10);
	    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(7,8,12);
	    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(5,5,5);
	    EXPECT_EQ(0, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(2,3,3);
	    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
	    vec.setAllCoordinates(30,17,13);
	    EXPECT_EQ(1, ingredients.getLatticeEntry(vec));
}