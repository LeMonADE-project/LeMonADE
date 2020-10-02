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

#include <LeMonADE/updater/moves/MoveConnectBase.h>


class TestMoveConnectSc: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_3( FeatureMoleculesIO, FeatureConnectionSc,FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<bool> >) Features;
  typedef ConfigureSystem<VectorInt3,Features,3> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;
  const IngredientsType& getIngredients() const {return ingredients;}

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


TEST_F(TestMoveConnectSc, initialiseSetterGetter)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(1,ingredients.getMolecules().size());

  MoveConnectSc move;
  
  move.init(ingredients);
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(2,move.getDir().getLength());
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(),move.getPartner());
  
  move.init(ingredients,0);
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(2,move.getDir().getLength());
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(),move.getPartner());
  
  
  VectorInt3 dir(2,0,0);
  move.init(ingredients,0,dir);
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(),move.getPartner());
  EXPECT_EQ(dir,move.getDir());
  
  VectorInt3 dir2(2,0,1);
  EXPECT_ANY_THROW(move.init(ingredients,0,dir2));
  ingredients.modifyMolecules().addMonomer(10,8,8);
  EXPECT_NO_THROW(ingredients.synchronize());
  move.init(ingredients,0,dir);
  EXPECT_EQ(1,move.getPartner());
}


TEST_F(TestMoveConnectSc, apply)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(6,8,8);
  ingredients.modifyMolecules().addMonomer(10,8,8);
  ingredients.modifyMolecules().addMonomer(8,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(3,ingredients.getMolecules().size());

  MoveConnectSc move;
  VectorInt3 dir(2,0,0);
  move.init(ingredients,2,dir);
  move.apply(ingredients); //first connection 
  
  dir.setAllCoordinates(-2,0,0);
  move.init(ingredients,2,dir);
  move.apply(ingredients); //second connection 
  
  EXPECT_TRUE(ingredients.getMolecules().areConnected(2,1));
  EXPECT_TRUE(ingredients.getMolecules().areConnected(2,0));
  EXPECT_FALSE(ingredients.getMolecules().areConnected(1,0));
  
}
