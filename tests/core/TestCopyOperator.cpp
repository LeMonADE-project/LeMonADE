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

/**
 * @file
 * @brief Tests for copy constructors
 *
 * @author Martin
 * @date 10.05.2015
 * */
/*****************************************************************************/

#include "gtest/gtest.h"
#include<stdint.h>

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>

using namespace std;

class CopyOperatorTest: public ::testing::Test{
  /* suppress cout output for better readability -->un-/comment here: */
public:
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
  /* ** */
};


TEST_F(CopyOperatorTest,Molecules)
{
    typedef LOKI_TYPELIST_1(FeatureMoleculesIO) Features1;
    typedef ConfigureSystem<VectorInt3,Features1,4> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 ingredients;

    //  #### prepare ingredients ####
    ingredients.setBoxX(64);
    ingredients.setBoxY(64);
    ingredients.setBoxZ(64);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(0);

    //bondset
    //ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.modifyBondset().addBond(2,1,0,23);
    ingredients.modifyBondset().addBond(-2,-1,0,32);
    ingredients.synchronize(ingredients);

    //molecules: size, positions, maxConnectivity and connections
    ingredients.modifyMolecules().resize(10);
    for(uint32_t i=0;i<10;i++){
      ingredients.modifyMolecules()[i].setAllCoordinates(2*i+10,i+10,10);
      if( (i>0) && (i<6) ){
	ingredients.modifyMolecules().connect(i-1,i,3);
      }
    }
    //compressed=solvent indizees
    //ingredients.setCompressedOutputIndices(6,9);
    //set age
    ingredients.modifyMolecules().setAge(42);

    std::cout<< "preparation of ingredients finished\n";
    ingredients.synchronize(ingredients);
    std::cout<< "synchronization of ingredients finished\n";


    Ing1 copyOfIngredients(ingredients);
    Ing1 assignedIngredients;
    assignedIngredients=ingredients;

    //check size of molecules
    EXPECT_EQ(ingredients.getMolecules().size(),copyOfIngredients.getMolecules().size());
    EXPECT_EQ(ingredients.getMolecules().size(),assignedIngredients.getMolecules().size());
    //check max connectivity
    EXPECT_EQ(ingredients.getMolecules().getMaxConnectivity(),copyOfIngredients.getMolecules().getMaxConnectivity());
    EXPECT_EQ(ingredients.getMolecules().getMaxConnectivity(),assignedIngredients.getMolecules().getMaxConnectivity());
    //check age
    EXPECT_EQ(ingredients.getMolecules().getAge(),copyOfIngredients.getMolecules().getAge());
    EXPECT_EQ(ingredients.getMolecules().getAge(),assignedIngredients.getMolecules().getAge());
    //check num Links
    EXPECT_EQ(ingredients.getMolecules().getTotalNumLinks(),copyOfIngredients.getMolecules().getTotalNumLinks());
    EXPECT_EQ(ingredients.getMolecules().getTotalNumLinks(),assignedIngredients.getMolecules().getTotalNumLinks());
    //check compressed indizees
    EXPECT_EQ(ingredients.getCompressedOutputIndices(),copyOfIngredients.getCompressedOutputIndices());
    EXPECT_EQ(ingredients.getCompressedOutputIndices(),assignedIngredients.getCompressedOutputIndices());
    //check positions
    for(uint i=0;i<ingredients.getMolecules().size();i++){
      EXPECT_EQ(ingredients.getMolecules()[i],copyOfIngredients.getMolecules()[i]);
      EXPECT_EQ(ingredients.getMolecules()[i],assignedIngredients.getMolecules()[i]);
    }
    //check connections
    for(uint i=0;i<ingredients.getMolecules().size();i++){
      for(uint j=0;j<ingredients.getMolecules().getNumLinks(i);j++){
	EXPECT_EQ(ingredients.getMolecules().getNeighborIdx(i,j),copyOfIngredients.getMolecules().getNeighborIdx(i,j));	EXPECT_EQ(ingredients.getMolecules().getNeighborIdx(i,j),assignedIngredients.getMolecules().getNeighborIdx(i,j));
      }
    }
    //check edges
    for(uint i=0;i<ingredients.getMolecules().size();i++){
      for(uint j=i;j<ingredients.getMolecules().size();j++){
	//check adges
	if(ingredients.getMolecules().areConnected(i,j)){
	  EXPECT_EQ(ingredients.getMolecules().getLinkInfo(i,j), copyOfIngredients.getMolecules().getLinkInfo(i,j));
	  EXPECT_EQ(ingredients.getMolecules().getLinkInfo(i,j), assignedIngredients.getMolecules().getLinkInfo(i,j));
	}else{
	  EXPECT_THROW(copyOfIngredients.getMolecules().getLinkInfo(i,j),std::runtime_error);
	  EXPECT_THROW(assignedIngredients.getMolecules().getLinkInfo(i,j),std::runtime_error);
	}
      }
    }

    //check box
    EXPECT_EQ(ingredients.getBoxX(),copyOfIngredients.getBoxX());
    EXPECT_EQ(ingredients.getBoxY(),copyOfIngredients.getBoxY());
    EXPECT_EQ(ingredients.getBoxZ(),copyOfIngredients.getBoxZ());
    EXPECT_EQ(ingredients.getBoxX(),assignedIngredients.getBoxX());
    EXPECT_EQ(ingredients.getBoxY(),assignedIngredients.getBoxY());
    EXPECT_EQ(ingredients.getBoxZ(),assignedIngredients.getBoxZ());
    //check periodicity
    EXPECT_EQ(ingredients.isPeriodicX(),copyOfIngredients.isPeriodicX());
    EXPECT_EQ(ingredients.isPeriodicY(),copyOfIngredients.isPeriodicY());
    EXPECT_EQ(ingredients.isPeriodicZ(),copyOfIngredients.isPeriodicZ());
    EXPECT_EQ(ingredients.isPeriodicX(),assignedIngredients.isPeriodicX());
    EXPECT_EQ(ingredients.isPeriodicY(),assignedIngredients.isPeriodicY());
    EXPECT_EQ(ingredients.isPeriodicZ(),assignedIngredients.isPeriodicZ());
    //check bondset
    EXPECT_EQ(ingredients.getBondset().getBondVector(23),copyOfIngredients.getBondset().getBondVector(23));
    EXPECT_EQ(ingredients.getBondset().getBondVector(23),assignedIngredients.getBondset().getBondVector(23));
    EXPECT_EQ(ingredients.getBondset().getBondVector(32),copyOfIngredients.getBondset().getBondVector(32));
    EXPECT_EQ(ingredients.getBondset().getBondVector(32),assignedIngredients.getBondset().getBondVector(32));

}


TEST_F(CopyOperatorTest,FeatureExcludedVolumeSc)
{
    typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureExcludedVolumeSc< FeatureLattice <int8_t> >) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 ingredients;

    //prepare ingredients
    ingredients.setBoxX(16);
    ingredients.setBoxY(16);
    ingredients.setBoxZ(16);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);

    //bondset
    //ingredients.modifyBondset().addBond(2,1,0,23);
    //ingredients.synchronize(ingredients);

    ingredients.modifyMolecules().resize(4);
    for(uint32_t i=0;i<4;i++){
      ingredients.modifyMolecules()[i].setAllCoordinates(2*i+1,i+1,1);
    }
    ingredients.synchronize(ingredients);

    Ing1 copyOfIngredients(ingredients);
    Ing1 assignedIngredients;
    assignedIngredients=ingredients;

    for(uint i=0;i<16;i++){
      for(uint j=0;j<16;j++){
	for(uint k=0;k<16;k++){
	  EXPECT_EQ(ingredients.getLatticeEntry(i,j,k),copyOfIngredients.getLatticeEntry(i,j,k));
	  EXPECT_EQ(ingredients.getLatticeEntry(i,j,k),assignedIngredients.getLatticeEntry(i,j,k));
	}
      }
    }
    MoveLocalSc testmove1,testmove2,testmove3;
    while((testmove1.getDir().getX()!=1) || (testmove1.getIndex()!=1)) testmove1.init(ingredients);
    while((testmove2.getDir().getX()!=1) || (testmove2.getIndex()!=1)) testmove2.init(copyOfIngredients);
    while((testmove3.getDir().getX()!=1) || (testmove3.getIndex()!=1)) testmove3.init(assignedIngredients);
    EXPECT_EQ(testmove1.check(ingredients),testmove2.check(copyOfIngredients));
    EXPECT_EQ(testmove1.check(ingredients),testmove3.check(assignedIngredients));
    EXPECT_EQ(testmove1.check(ingredients),testmove1.check(copyOfIngredients));
    EXPECT_EQ(testmove1.check(ingredients),testmove1.check(assignedIngredients));
}


TEST_F(CopyOperatorTest,FeatureAttributes)
{
    typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureAttributes<>) Features1;
    typedef ConfigureSystem<VectorInt3,Features1> Config1;
    typedef Ingredients<Config1> Ing1;
    Ing1 ingredients;

    //prepare ingredients
    ingredients.setBoxX(16);
    ingredients.setBoxY(16);
    ingredients.setBoxZ(16);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);

    ingredients.modifyMolecules().resize(4);
    for(uint32_t i=0;i<4;i++){
      ingredients.modifyMolecules()[i].setAllCoordinates(2*i+1,i+1,1);
      ingredients.modifyMolecules()[i].setAttributeTag(i+2*i);
    }
    ingredients.synchronize(ingredients);

    Ing1 copyOfIngredients(ingredients);
    Ing1 assignedIngredients;
    assignedIngredients=ingredients;

    for(uint32_t i=0;i<4;i++){
      EXPECT_EQ(ingredients.getMolecules()[i].getAttributeTag(),copyOfIngredients.getMolecules()[i].getAttributeTag());
      EXPECT_EQ(ingredients.getMolecules()[i].getAttributeTag(),assignedIngredients.getMolecules()[i].getAttributeTag());
    }

}
