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
 * @brief Tests for the class MoveLocalSc
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>

#include <LeMonADE/updater/moves/MoveLocalSc.h>


class TestMoveLocalSc: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<bool> >) Features;
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

TEST_F(TestMoveLocalSc, initialiseSetterGetter)
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

  MoveLocalSc move;

  // ######################################################################## //
  // use empty init interface: dice a random monomer and a random move direction
  move.init(getIngredients());
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(1,move.getDir().getLength());
  //change probability
  move.multiplyProbability(0.5);

  //add a new monomer, check if init changes properties correctly
  ingredients.modifyMolecules().addMonomer(4,8,8);
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(2,ingredients.getMolecules().size());

  move.init(getIngredients());
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(1,move.getDir().getLength());

  // ######################################################################## //
  // use index init interface: set the index and dice a random move direction

  EXPECT_NO_THROW(move.init(getIngredients(),0));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(1,move.getDir().getLength());
  //change probability
  move.multiplyProbability(0.5);

  EXPECT_EQ(2,ingredients.getMolecules().size());
  EXPECT_NO_THROW(move.init(getIngredients(),1));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(1,move.getDir().getLength());

  //check wrong index by initialize
  EXPECT_ANY_THROW(move.init(getIngredients(),ingredients.getMolecules().size()));
  EXPECT_ANY_THROW(move.init(getIngredients(),2)); // ( =same as molecules.size() )
  EXPECT_ANY_THROW(move.init(getIngredients(),-1));
  EXPECT_ANY_THROW(move.init(getIngredients(),68468468));

  // ######################################################################## //
  // use direction init interface: dice an index and set the move direction
  VectorInt3 direction(1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());
  //change probability
  move.multiplyProbability(0.5);

  direction.setAllCoordinates(0,1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(0,0,1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(0,-1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(0,0,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_TRUE((0==move.getIndex()) || (1==move.getIndex()));
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(2,0,0);
  EXPECT_ANY_THROW(move.init(getIngredients(),direction));
  direction.setAllCoordinates(0,-2,0);
  EXPECT_ANY_THROW(move.init(getIngredients(),direction));

  // ######################################################################## //
  // use direction and index init interface: set index and the move direction
  //check all possible directions
  direction.setAllCoordinates(1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(0,move.getIndex());
  EXPECT_EQ(direction,move.getDir());
  //change probability
  move.multiplyProbability(0.5);

  direction.setAllCoordinates(0,1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(0,0,1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(-1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(0,-1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  direction.setAllCoordinates(0,0,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_EQ(1.0,move.getProbability());
  EXPECT_EQ(1,move.getIndex());
  EXPECT_EQ(direction,move.getDir());

  //check some forbidden directions and idicees
  EXPECT_ANY_THROW(move.init(getIngredients(),2,direction));
  EXPECT_ANY_THROW(move.init(getIngredients(),-1,direction));
  direction.setAllCoordinates(2,0,0);
  EXPECT_ANY_THROW(move.init(getIngredients(),1,direction));
  direction.setAllCoordinates(0,-1,1);
  EXPECT_ANY_THROW(move.init(getIngredients(),1,direction));
  EXPECT_ANY_THROW(move.init(getIngredients(),2,direction));

}

TEST_F(TestMoveLocalSc, checkAndApply)
{
  // starting here with a new instance of ingredients!
  // class member lifetime is only one testcase!
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(false);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(2,8,8);
  ingredients.modifyMolecules().addMonomer(0,8,8);

  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(2,ingredients.getMolecules().size());

  MoveLocalSc move;

  // ################################# //
  //check move with FeatureExcludedVolumeSc:
  VectorInt3 direction(1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(0,1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(0,-1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(0,0,1);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(0,0,-1);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_TRUE(move.check(ingredients));

  direction.setAllCoordinates(1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients);
  // now monomers are at positions 3,8,8 and 0,8,8

  // ################################# //
  //check move with FeatureMoleculesIO:
  //bondset:
  ingredients.modifyMolecules().connect(0,1);
  direction.setAllCoordinates(1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(-1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),0,direction));
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients);
  EXPECT_EQ(VectorInt3(2,8,8),VectorInt3(ingredients.getMolecules()[0]));
  EXPECT_EQ(VectorInt3(0,8,8),VectorInt3(ingredients.getMolecules()[1]));
  // now monomers are at positions 2,8,8 and 0,8,8

  //nonperiodic wall
  direction.setAllCoordinates(-1,0,0);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_FALSE(move.check(ingredients));

  // ################################# //
  //check move with FeatureFixedMonomers:
  ingredients.modifyMolecules()[0].setMovableTag(false);
  EXPECT_NO_THROW(move.init(getIngredients(),0));
  EXPECT_FALSE(move.check(ingredients));

  direction.setAllCoordinates(0,1,0);
  EXPECT_NO_THROW(move.init(getIngredients(),1,direction));
  EXPECT_TRUE(move.check(ingredients));

  move.apply(ingredients);

  EXPECT_EQ(VectorInt3(2,8,8),VectorInt3(ingredients.getMolecules()[0]));
  EXPECT_EQ(VectorInt3(0,9,8),VectorInt3(ingredients.getMolecules()[1]));

}
