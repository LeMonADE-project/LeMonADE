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
 * @brief Tests for the class MoveAddMonomerSc
 * @author Martin
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureAttributes.h>

#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>


class TestMoveAddMonomerSc: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes<>, FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<bool> >) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
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

TEST_F(TestMoveAddMonomerSc, initialiseSetterGetter)
{
  //IngredientsType ingredients;
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());

  MoveAddMonomerSc<> addmove;

  // init: set particle index to size()+1, set probability to 0
  addmove.init(ingredients);

  //check particle index
  EXPECT_EQ(1,addmove.getMonomerIndex());
  //check probability
  EXPECT_DOUBLE_EQ(1.0,addmove.getProbability());
  //change and reset probability
  addmove.multiplyProbability(0.5);
  EXPECT_DOUBLE_EQ(0.5,addmove.getProbability());
  addmove.resetProbability();
  EXPECT_DOUBLE_EQ(1.0,addmove.getProbability());

  //check setter/getter for type
  addmove.setTag(3);
  EXPECT_EQ(3,addmove.getTag());

  //check overwriting of type
  addmove.setTag(9);
  EXPECT_EQ(9,addmove.getTag());

  //check setter/getter of position VectorInt3
  VectorInt3 position(0,8,8);
  addmove.setPosition(position);
  EXPECT_EQ(position, addmove.getPosition());

  //check overwriting of position VectorInt3
  position.setAllCoordinates(8,8,0);
  addmove.setPosition(position);
  EXPECT_EQ(position, addmove.getPosition());

  //check setter/getter of position (int, int, int)
  addmove.setPosition(8,0,8);
  EXPECT_EQ(VectorInt3(8,0,8),addmove.getPosition());

  //check overwriting of position (int, int, int)
  addmove.setPosition(0,0,8);
  EXPECT_EQ(VectorInt3(0,0,8),addmove.getPosition());

}

TEST_F(TestMoveAddMonomerSc, checkAndApply)
{
  //IngredientsType ingredients;
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(false);
  ingredients.setPeriodicY(false);
  ingredients.setPeriodicZ(false);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());

  MoveAddMonomerSc<> addmove;

  // init: set particle index to size()+1, set probability to 0
  addmove.init(ingredients);
  addmove.setTag(5);
  addmove.setPosition(8,8,8);

  //check FeatureExcludedVolumeSc (all lattice positions of monomer 0)
  addmove.setPosition(8,8,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(9,8,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(8,9,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(8,8,9);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(9,9,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(9,8,9);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(8,9,9);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(9,9,9);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(7,8,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(8,7,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(8,8,7);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(7,7,8);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(8,7,7);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(7,7,7);
  EXPECT_FALSE(addmove.check(ingredients));

  //check some free positions
  addmove.setPosition(8,8,6);
  EXPECT_TRUE(addmove.check(ingredients));
  addmove.setPosition(8,7,6);
  EXPECT_TRUE(addmove.check(ingredients));
  addmove.setPosition(6,8,7);
  EXPECT_TRUE(addmove.check(ingredients));
  addmove.setPosition(7,8,6);
  EXPECT_TRUE(addmove.check(ingredients));
  addmove.setPosition(9,6,8);
  EXPECT_TRUE(addmove.check(ingredients));

  //check box boundaries
  addmove.setPosition(0,0,15);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(0,15,0);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(15,0,0);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(0,15,15);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(15,15,0);
  EXPECT_FALSE(addmove.check(ingredients));
  addmove.setPosition(15,15,15);
  EXPECT_FALSE(addmove.check(ingredients));

  // check corner position
  addmove.setPosition(0,0,0);
  EXPECT_TRUE(addmove.check(ingredients));

  //apply a move and check consistency
  addmove.apply(ingredients);
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(2,ingredients.getMolecules().size());
  EXPECT_EQ(VectorInt3(0,0,0),ingredients.getMolecules()[1]);
  EXPECT_EQ(5,ingredients.getMolecules()[1].getAttributeTag());
}
