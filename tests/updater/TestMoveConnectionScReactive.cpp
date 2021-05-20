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

/*****************************************************************************/
/**
 * @file
 * @brief Tests for the class MoveConnectSc
 * @author Toni 
 * */
/*****************************************************************************/

#include <limits>

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>

#include <LeMonADE/updater/moves/MoveConnectScReactive.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>


class TestMoveConnectionScReactive: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_3( FeatureMoleculesIO,FeatureReactiveBonds, FeatureExcludedVolumeSc<>) Features;
  typedef ConfigureSystem<VectorInt3,Features,3> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;

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


TEST_F(TestMoveConnectionScReactive, checkAll)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);
  ingredients.modifyMolecules().addMonomer(10,8,8);
  ingredients.modifyMolecules().addMonomer(10,10,8);
  ingredients.modifyMolecules()[0].setReactive(true);
  ingredients.modifyMolecules()[1].setReactive(true);
  ingredients.modifyMolecules()[2].setReactive(false);
  ingredients.modifyMolecules()[0].setNumMaxLinks(1);
  ingredients.modifyMolecules()[1].setNumMaxLinks(2);
  ingredients.modifyMolecules()[2].setNumMaxLinks(2);
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(3,ingredients.getMolecules().size());
  EXPECT_FALSE(ingredients.getMolecules().areConnected(0,1));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(1,2));
  
  MoveConnectScReactive move;
  
  move.init(ingredients);
  EXPECT_TRUE((move.getIndex() >= 0) && (move.getIndex() <= 1)); // only 0 and 1 are reactive
  EXPECT_TRUE((move.getDir()*move.getDir() >= 4) && (move.getDir()*move.getDir() <= 6));
  // will maybe fail as the direction is selected randomly hence the partner also
  //EXPECT_TRUE((move.getPartner() >= 0) && (move.getPartner() <= 1)); // only 0 and 1 are reactive
  EXPECT_TRUE(((move.getPartner() >= 0) && (move.getPartner() <= 1)) || (move.getPartner()==std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  EXPECT_FALSE((move.getPartner() >= 2) && (move.getPartner() < std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  
  
  move.init(ingredients,0);
  EXPECT_EQ(0, move.getIndex());
  EXPECT_TRUE((move.getIndex() >= 0) && (move.getIndex() <= 1)); // only 0 and 1 are reactive
  EXPECT_TRUE((move.getDir()*move.getDir() >= 4) && (move.getDir()*move.getDir() <= 6));
  EXPECT_TRUE(((move.getPartner() == 0)) || (move.getPartner()==std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  EXPECT_FALSE((move.getPartner() >= 2) && (move.getPartner() < std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  

  move.init(ingredients,1);
  EXPECT_EQ(1, move.getIndex());
  EXPECT_TRUE((move.getIndex() >= 0) && (move.getIndex() <= 1)); // only 0 and 1 are reactive
  EXPECT_TRUE((move.getDir()*move.getDir() >= 4) && (move.getDir()*move.getDir() <= 6));
  EXPECT_TRUE(((move.getPartner() == 0)) || (move.getPartner()==std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  EXPECT_FALSE((move.getPartner() >= 2) && (move.getPartner() < std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  
  
  move.init(ingredients,1,2);
  EXPECT_EQ(move.getDir().getY(),2);
  EXPECT_FALSE(move.check(ingredients));
  EXPECT_TRUE(move.getIndex() == 1); // only 0 and 1 are reactive
  EXPECT_TRUE((move.getDir()*move.getDir() == 4));
  EXPECT_TRUE((move.getPartner() == 2)); // only 0 and 1 are reactive
  EXPECT_FALSE((move.getPartner() < 2) && (move.getPartner() < std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  
  
  move.init(ingredients,0,1);
  EXPECT_EQ(move.getDir().getX(),2);
  EXPECT_TRUE(move.getIndex() == 0); // only 0 and 1 are reactive
  EXPECT_TRUE((move.getDir()*move.getDir() == 4));
  EXPECT_TRUE((move.getPartner() == 1)); // only 0 and 1 are reactive
  EXPECT_FALSE((move.getPartner() < 1) && (move.getPartner() < std::numeric_limits<uint32_t>::max()) ); // only 0 and 1 are reactive
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients);
  
  EXPECT_TRUE(ingredients.getMolecules().areConnected(0,1));
  
  move.init(ingredients,0,1); //already occupied
  EXPECT_FALSE(move.check(ingredients));
  
  //no reactive groups left
  EXPECT_EQ(ingredients.getNUnreactedMonomers(), 0);
  move.init(ingredients);
  EXPECT_EQ(move.getIndex(), 0 ); // default value for invalid move
  EXPECT_TRUE(move.getPartner() == std::numeric_limits<uint32_t>::max());
  EXPECT_FALSE(move.check(ingredients));
  

 
}
