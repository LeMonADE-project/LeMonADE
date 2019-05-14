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

/************************************************************************/
/**
 * @file
 * @brief Tests for the bondset feature
 *
 * Includes tests for the following classes:  SlowBondset
 * */
/************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/utility/SlowBondset.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/GenerateMonomerType.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class SlowBondsetTest: public ::testing::Test{
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

/************************************************************************/
//test the addBond methods of class SlowBondset
/************************************************************************/
TEST_F(SlowBondsetTest, AddBonds){
  SlowBondset bondset;

  EXPECT_EQ(0,bondset.size());
  //add some vectors
  bondset.addBond(2,2,1,77);

  EXPECT_EQ(1,bondset.size());

  VectorInt3 bondvector(2,0,0);
  int32_t identifier=17;

  bondset.addBond(bondvector, identifier);

  EXPECT_EQ(2,bondset.size());
  //try to add duplicate bondvector (same components)
  EXPECT_THROW(bondset.addBond(2,2,1,37),std::runtime_error);

  //try to add duplicate bondvector (same identifier)
  EXPECT_THROW(bondset.addBond(3,1,0,17),std::runtime_error);

}


/************************************************************************/
//test the copy constructor of class SlowBondset
/************************************************************************/
TEST_F(SlowBondsetTest, CopyConstructor){
  SlowBondset bondset;

  EXPECT_EQ(0,bondset.size());
  //add some vectors
  bondset.addBond(2,2,1,77);
  bondset.addBond(2,1,2,78);

  EXPECT_EQ(2,bondset.size());

  //copy bondset
  SlowBondset bondset_new(bondset);

  //check size
  EXPECT_EQ(2,bondset_new.size());

  //try to add duplicate bondvector (same components)
  EXPECT_THROW(bondset_new.addBond(2,2,1,37),std::runtime_error);

  //try to add duplicate bondvector (same identifier)
  EXPECT_THROW(bondset_new.addBond(3,1,0,78),std::runtime_error);

  //check lookup table
  VectorInt3 bondvector77(2,2,1);
  EXPECT_TRUE(bondset_new.isValid(bondvector77));
  EXPECT_TRUE(bondset_new.isValidStrongCheck(bondvector77));

  //check components
  EXPECT_EQ(bondset_new.getBondVector(77).getX(),2);
  EXPECT_EQ(bondset_new.getBondVector(77).getY(),2);
  EXPECT_EQ(bondset_new.getBondVector(77).getZ(),1);

  //check if not only referenced to old bondset
  bondset.addBond(1,2,2,79);
  VectorInt3 bondvector79(1,2,2);
  EXPECT_FALSE(bondset_new.isValid(bondvector79));
  EXPECT_FALSE(bondset_new.isValidStrongCheck(bondvector79));
  bondset.clear();
  EXPECT_EQ(0,bondset.size());
  EXPECT_EQ(2,bondset_new.size());
  EXPECT_TRUE(bondset_new.isValid(bondvector77));
  EXPECT_TRUE(bondset_new.isValidStrongCheck(bondvector77));
  EXPECT_EQ(bondset_new.getBondVector(77).getX(),2);
  EXPECT_EQ(bondset_new.getBondVector(77).getY(),2);
  EXPECT_EQ(bondset_new.getBondVector(77).getZ(),1);

}

/************************************************************************/
//test the getBondIdentifier and getBondVector methods of class SlowBondset
/************************************************************************/
TEST_F(SlowBondsetTest,GetBondInformation){
  SlowBondset bondset;

  //add some vectors
  bondset.addBond(2,2,1,77);
  bondset.addBond(-2,0,0,20);


  //check reaction to invalid requests for vectors and identifiers
  EXPECT_THROW(bondset.getBondIdentifier(2,2,-1),std::runtime_error);
  EXPECT_THROW(bondset.getBondVector(30),std::runtime_error);

  //check if getter return the correct values

  //check returned identifiers
  EXPECT_EQ(bondset.getBondIdentifier(2,2,1),77);
  EXPECT_EQ(bondset.getBondIdentifier(-2,0,0),20);

  //check returned vectors
  //metadata component
  //EXPECT_EQ(bondset.getBondVector(77).getW(),77);
  //EXPECT_EQ(bondset.getBondVector(20).getW(),20);
  //spacial components
  EXPECT_EQ(bondset.getBondVector(20).getX(),-2);
  EXPECT_EQ(bondset.getBondVector(20).getY(),0);
  EXPECT_EQ(bondset.getBondVector(20).getZ(),0);

  EXPECT_EQ(bondset.getBondVector(77).getX(),2);
  EXPECT_EQ(bondset.getBondVector(77).getY(),2);
  EXPECT_EQ(bondset.getBondVector(77).getZ(),1);
}

TEST_F(SlowBondsetTest, LookupTable)
{

  SlowBondset bondset;

  //set up some vectors
  VectorInt3 v1(1,0,-3);
  VectorInt3 v2(-1,2,2);
  VectorInt3 v3(2,2,1);
  VectorInt3 v4(1,1,1);
  VectorInt3 v5(0,0,0);
  VectorInt3 v6(4,0,0);
  VectorInt3 v7(1,-3,0);

  //checking bonds before lookup was synchronized should cause an exception
  EXPECT_THROW(bondset.isValid(v1),std::runtime_error);

  //resetting the bondset should not change this behaviour
  EXPECT_THROW(bondset.isValid(v2),std::runtime_error);

  //updating the lookup now without any bondvectors set should cause all vectors to be invalid
  bondset.updateLookupTable();

  EXPECT_FALSE(bondset.isValid(v1));
  EXPECT_FALSE(bondset.isValid(v2));
  EXPECT_FALSE(bondset.isValid(v3));

  EXPECT_FALSE(bondset.isValid(v4));

  EXPECT_FALSE(bondset.isValid(v5));

  EXPECT_FALSE(bondset.isValid(v6));

  EXPECT_FALSE(bondset.isValid(v7));

  //resetting the lookup should make the previous bahaviour return
  bondset.resetLookupTable();
  EXPECT_THROW(bondset.isValid(v2),std::runtime_error);

  //adding a bondvector but not updating the lookup should cause an exception when checking a vector
  bondset.addBond(v1,105);
  bondset.addBond(-1,2,2,91);
  bondset.addBond(v3,77);
  EXPECT_THROW(bondset.isValid(v1),std::runtime_error);

  bondset.updateLookupTable();
  //now the added bonds should be valid, the other ones not
  EXPECT_TRUE(bondset.isValid(v1));
  EXPECT_TRUE(bondset.isValid(v2));
  EXPECT_TRUE(bondset.isValid(v3));

  EXPECT_FALSE(bondset.isValid(v4));
  EXPECT_FALSE(bondset.isValid(v5));
  EXPECT_FALSE(bondset.isValid(v6));
  EXPECT_FALSE(bondset.isValid(v7));

  //again resetting the lookup should bring back the old behaviour
  bondset.resetLookupTable();
  EXPECT_THROW(bondset.isValid(v2),std::runtime_error);

}
