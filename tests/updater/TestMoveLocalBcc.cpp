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
 * @brief Tests for the class MoveLocalBcc
 * @author Martin
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeBcc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>

#include <LeMonADE/updater/moves/MoveLocalBcc.h>


class TestMoveLocalBcc: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureExcludedVolumeBcc<FeatureLatticePowerOfTwo<bool> >) Features;
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

TEST_F(TestMoveLocalBcc, initialiseSetterGetter)
{
  //IngredientsType ingredients;
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBccBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());

  MoveLocalBcc move;

  // ######################################################################## //
  // use empty init interface: dice a random monomer and a random move direction
  move.init(getIngredients());
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_DOUBLE_EQ(sqrt(3),move.getDir().getLength());
  //change probability
  move.multiplyProbability(0.5);

  //add a new monomer, check if init changes properties correctly
  ingredients.modifyMolecules().addMonomer(4,8,8);
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(2,ingredients.getMolecules().size());

  move.init(getIngredients());
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_DOUBLE_EQ(sqrt(3),move.getDir().getLength());

  // ######################################################################## //
  // use index init interface: dice a direction and set the index

  //check correct index by initialize
  EXPECT_NO_THROW(move.init(getIngredients(),0));
  EXPECT_EQ(0,move.getIndex());
  EXPECT_NO_THROW(move.init(getIngredients(),1));
  EXPECT_EQ(1,move.getIndex());

  //check wrong index by initialize
  EXPECT_ANY_THROW(move.init(getIngredients(),ingredients.getMolecules().size()));
  EXPECT_ANY_THROW(move.init(getIngredients(),2)); // ( =same as molecules.size() )
  EXPECT_ANY_THROW(move.init(getIngredients(),-1));
  EXPECT_ANY_THROW(move.init(getIngredients(),68468468));

  // ######################################################################## //
  // use direction init interface: dice an index and set the move direction
  VectorInt3 direction(1,1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());
  //change probability
  move.multiplyProbability(0.5);

  direction.setAllCoordinates(1,1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(1,-1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,-1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(1,-1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,-1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(2,0,0);
  EXPECT_ANY_THROW(move.init(getIngredients(),direction));
  direction.setAllCoordinates(0,-1,0);
  EXPECT_ANY_THROW(move.init(getIngredients(),direction));

  // ######################################################################## //
  // use direction and index init interface: set index and the move direction
  direction.setAllCoordinates(1,1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(direction,move.getDir());
  //change probability
  move.multiplyProbability(0.5);

  direction.setAllCoordinates(1,1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(1,-1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(1,-1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,-1,1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,-1,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  EXPECT_ANY_THROW(move.init(getIngredients(),2,direction));
  EXPECT_ANY_THROW(move.init(getIngredients(),-1,direction));
  direction.setAllCoordinates(2,1,-1);
  EXPECT_ANY_THROW(move.init(getIngredients(),1,direction));
  direction.setAllCoordinates(0,-1,1);
  EXPECT_ANY_THROW(move.init(getIngredients(),1,direction));
  EXPECT_ANY_THROW(move.init(getIngredients(),2,direction));
}

TEST_F(TestMoveLocalBcc, checkAndApply)
{
  // starting here with a new instance of ingredients!
  // class member lifetime is only one testcase!
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(false);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBccBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(8,8,8);
  ingredients.modifyMolecules().addMonomer(6,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());

  MoveLocalBcc move;

  // ######################################################################## //
  // check feature Excluded Volume
  VectorInt3 direction(1,1,1);
  move.init(ingredients,1,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(1,1,-1);
  move.init(ingredients,1,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(1,-1,1);
  move.init(ingredients,1,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(1,-1,-1);
  move.init(ingredients,1,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,-1);
  move.init(ingredients,0,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,1,-1);
  move.init(ingredients,0,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,1);
  move.init(ingredients,0,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,1,1);
  move.init(ingredients,0,direction);
  EXPECT_FALSE(move.check(ingredients));

  //check some possible moves
  direction.setAllCoordinates(-1,1,1);
  move.init(ingredients,1,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,1);
  move.init(ingredients,1,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(-1,1,-1);
  move.init(ingredients,1,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,-1);
  move.init(ingredients,1,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,1,1);
  move.init(ingredients,0,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,-1,1);
  move.init(ingredients,0,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,1,-1);
  move.init(ingredients,0,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,-1,-1);
  move.init(ingredients,0,direction);
  EXPECT_TRUE(move.check(ingredients));

  move.apply(ingredients);
  EXPECT_NO_THROW(ingredients.synchronize());

  // ######################################################################## //
  // check feature Box
  ingredients.modifyMolecules().addMonomer(0,0,0);
  EXPECT_NO_THROW(ingredients.synchronize());

  direction.setAllCoordinates(-1,1,1);
  move.init(ingredients,2,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,1);
  move.init(ingredients,2,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,1,-1);
  move.init(ingredients,2,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,-1,-1);
  move.init(ingredients,2,direction);
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(1,1,1);
  move.init(ingredients,2,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,-1,1);
  move.init(ingredients,2,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,1,-1);
  move.init(ingredients,2,direction);
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,-1,-1);
  move.init(ingredients,2,direction);
  EXPECT_TRUE(move.check(ingredients));

  //apply move and check again
  move.apply(ingredients);
  EXPECT_NO_THROW(ingredients.synchronize());
  direction.setAllCoordinates(-1,1,1);
  move.init(ingredients,2,direction);
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients);
  direction.setAllCoordinates(-1,-1,-1);
  move.init(ingredients,2,direction);
  EXPECT_FALSE(move.check(ingredients));
}
