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

/*****************************************************************************/
/**
 * @file
 * @brief Tests for UpdaterAbstractCreate
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/updater/UpdaterAbstractCreate.h>

using namespace std;




/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class UpdaterAbstractCreateTest: public ::testing::Test{
public:

  //redirect cout output
  virtual void SetUp(){
    originalBuffer=cout.rdbuf();
    cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    cout.rdbuf(originalBuffer);
  };

private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;



};

//Define a test class of the updater
template<class IngredientsType>
class UpdaterTestAbstractCreate: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
  UpdaterTestAbstractCreate(IngredientsType& ingredients_):BaseClass(ingredients_),numExec(0){}

  virtual void initialize(){}
  virtual bool execute(){
    bool checker(true);
    if(numExec==0){
      checker=addSingleMonomer(1);
      checker=addSingleMonomer(2);
      checker=addSingleMonomer(3);
    }else if(numExec==1){
      ingredients.modifyMolecules().resize(0);
      ingredients.synchronize();

      checker=addMonomerAtPosition(VectorInt3( ingredients.getBoxX()/2,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),4);
      checker=addMonomerToParent(0,3);
    }else if(numExec==2){
      checker=addMonomerToParent(1,3);
    }else if(numExec==3){
      moveSystem(20);
    }else if(numExec==4){
      checker=addMonomerInsideConnectedPair(1,2,4);
      checker=addMonomerInsideConnectedPair(1,ingredients.getMolecules().size()-1,4);
      checker=addMonomerInsideConnectedPair(1,ingredients.getMolecules().size()-1,4);
    }else if(numExec==5){
      linearizeSystem();
    }else if(numExec==6){
      ingredients.modifyMolecules().resize(0);
      ingredients.synchronize();
      //setup test chain
      checker=addMonomerAtPosition(VectorInt3(2,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),1);
      checker=addMonomerAtPosition(VectorInt3(4,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),1);
      ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
      checker=addMonomerAtPosition(VectorInt3(6,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),1);
      ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
      checker=addMonomerAtPosition(VectorInt3(8,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),1);
      ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
      checker=addMonomerAtPosition(VectorInt3(10,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),1);
      ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
    }else if(numExec==7){
      //try adding a ring on non existing parent
      numExec++;
      checker=addRing(5);
    }else if(numExec==8){
      //add ring
      checker=addRing(2);
    }else if(numExec==9){
      //add ring
      checker=addRing(3,4,8);
    }
    numExec++;
    return checker;

  }
  virtual void cleanup(){}

  int32_t getNumExec(){return numExec;}

private:
  int32_t numExec;
  using BaseClass::ingredients;

  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::addMonomerAtPosition;
  using BaseClass::addMonomerInsideConnectedPair;
  using BaseClass::addRing;
  using BaseClass::moveSystem;
  using BaseClass::linearizeSystem;
};

TEST_F(UpdaterAbstractCreateTest, Constructor)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;

  //Constructor call
  EXPECT_NO_THROW(UpdaterAbstractCreate<IngredientsType> Timmy(ingredients));
}

TEST_F(UpdaterAbstractCreateTest, TestUpdater)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  EXPECT_NO_THROW(ingredients.synchronize());

  UpdaterTestAbstractCreate<IngredientsType> Tommy(ingredients);

  EXPECT_EQ(0,Tommy.getNumExec());
  //first execution
  EXPECT_TRUE(Tommy.execute());
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(3,ingredients.getMolecules().size());
  EXPECT_EQ(1,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[2].getAttributeTag());

  //second execution
  EXPECT_EQ(1,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(2,ingredients.getMolecules().size());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(8,ingredients.getMolecules()[0].getX());
  EXPECT_EQ(8,ingredients.getMolecules()[0].getY());
  EXPECT_EQ(8,ingredients.getMolecules()[0].getZ());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,0));
  EXPECT_EQ(2,(ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength());

  //second execution
  EXPECT_EQ(2,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(3,ingredients.getMolecules().size());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,2));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(2,1));
  EXPECT_EQ(2,(ingredients.getMolecules()[1]-ingredients.getMolecules()[2]).getLength());
  EXPECT_GE(3.17,(ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength());
  EXPECT_GE(5.17,(ingredients.getMolecules()[0]-ingredients.getMolecules()[2]).getLength());

  //third execution
  EXPECT_EQ(3,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(3,ingredients.getMolecules().size());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,2));
  EXPECT_GE(3.17,(ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength());
  EXPECT_GE(6.34,(ingredients.getMolecules()[0]-ingredients.getMolecules()[2]).getLength());

  //fourth execution
  EXPECT_EQ(4,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(6,ingredients.getMolecules().size());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[3].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[4].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[5].getAttributeTag());

  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,5));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(5,4));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(4,3));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(3,2));

  EXPECT_FALSE(ingredients.getMolecules().areConnected(0,5));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(0,4));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(0,3));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(0,2));

  EXPECT_FALSE(ingredients.getMolecules().areConnected(1,4));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(1,3));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(1,2));

  EXPECT_FALSE(ingredients.getMolecules().areConnected(5,0));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(5,3));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(5,2));

  EXPECT_FALSE(ingredients.getMolecules().areConnected(4,0));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(4,1));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(4,2));

  EXPECT_FALSE(ingredients.getMolecules().areConnected(3,0));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(3,1));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(3,5));

  EXPECT_FALSE(ingredients.getMolecules().areConnected(2,1));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(2,5));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(2,4));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(2,0));

  //fifth execution
  EXPECT_EQ(5,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(6,ingredients.getMolecules().size());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[3].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[4].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[5].getAttributeTag());

  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,2));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(2,3));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(3,4));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(4,5));

  //sixth execution
  EXPECT_EQ(6,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(5,ingredients.getMolecules().size());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,2));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(2,3));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(3,4));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(4,5));

  //seventh execution
  EXPECT_EQ(7, Tommy.getNumExec());
  EXPECT_ANY_THROW(Tommy.execute());
  EXPECT_EQ(5,ingredients.getMolecules().size());

  //Eigth execution 
  EXPECT_EQ(8, Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(10,ingredients.getMolecules().size());
  for (uint32_t j=0;j<5;j++)
  {
    EXPECT_EQ(1,ingredients.getMolecules()[5+j].getAttributeTag());
  }

  //Ninth execution 
  EXPECT_EQ(9, Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  EXPECT_EQ(18,ingredients.getMolecules().size());
  for (uint32_t j=0;j<8;j++)
  {
    EXPECT_EQ(4,ingredients.getMolecules()[10+j].getAttributeTag());
  }

}

