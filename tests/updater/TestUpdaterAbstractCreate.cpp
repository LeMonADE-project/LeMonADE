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

#include "UpdaterAbstractCreate.h"

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
      checker=add_lonely_monomer(1);
      checker=add_lonely_monomer(2);
      checker=add_lonely_monomer(3);
    }else if(numExec==1){
      ingredients.modifyMolecules().resize(0);
      ingredients.synchronize();
      
      checker=add_monomer_to_position(VectorInt3( ingredients.getBoxX()/2,ingredients.getBoxY()/2,ingredients.getBoxZ()/2),4);
      checker=add_monomer_to_parent(0,3);
    }else if(numExec==2){
      
    }
    
    numExec++;
    return checker;
    
  }
  virtual void cleanup(){}
  
  int32_t getNumExec(){return numExec;}
  
private:
  int32_t numExec;
  using BaseClass::ingredients;
  
  using BaseClass::rng;
  
  using BaseClass::add_monomer_to_parent;
  using BaseClass::add_lonely_monomer;
  using BaseClass::add_monomer_to_position;
  using BaseClass::move_system;
};

TEST_F(UpdaterAbstractCreateTest, Constructor)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;
  
  IngredientsType ingredients;
  
  //run the public functions
  UpdaterAbstractCreate<IngredientsType> Timmy(ingredients);
  
  //nothing should happen here
  Timmy.initialize();
  Timmy.execute();
  Timmy.cleanup();
}

TEST_F(UpdaterAbstractCreateTest, TestUpdater)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;
  
  IngredientsType ingredients;
  ingredients.setBoxX(8);
  ingredients.setBoxY(8);
  ingredients.setBoxZ(8);
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
  EXPECT_EQ(4,ingredients.getMolecules()[0].getX());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getY());
  EXPECT_EQ(4,ingredients.getMolecules()[0].getZ());
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(1,0));
  EXPECT_EQ(2,(ingredients.getMolecules()[0]-ingredients.getMolecules()[1]).getLength());
  
  //second execution
  EXPECT_EQ(2,Tommy.getNumExec());
  EXPECT_TRUE(Tommy.execute());
  
}

