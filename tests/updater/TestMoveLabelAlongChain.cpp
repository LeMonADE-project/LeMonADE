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
#include <LeMonADE/feature/FeatureLabel.h>

#include <LeMonADE/updater/moves/MoveLabelAlongChain.h>


class TestMoveLabelAlongChain: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureLabel) Features;
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

TEST_F(TestMoveLabelAlongChain, initialiseSetterGetter)
{
  //IngredientsType ingredients;
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(2,2,2);
  EXPECT_NO_THROW(ingredients.synchronize());

  
  MoveLabelAlongChain labelmove;

  // init: set particle index to size()+1, set probability to 0
  labelmove.init(ingredients);

  //check particle index
  EXPECT_EQ(0,labelmove.getIndex());
  //check probability
  EXPECT_DOUBLE_EQ(1.0,labelmove.getProbability());
  //change and reset probability
  labelmove.multiplyProbability(0.5);
  EXPECT_DOUBLE_EQ(0.5,labelmove.getProbability());
  labelmove.resetProbability();
  EXPECT_DOUBLE_EQ(1.0,labelmove.getProbability());

  //check setter/getter for dir
  labelmove.init(ingredients,0, MoveLabelAlongChain::RIGHT);
  EXPECT_EQ(1,labelmove.getDir());

  labelmove.init(ingredients,0, MoveLabelAlongChain::LEFT);
  EXPECT_EQ(-1,labelmove.getDir());

  //check setter and getter for label partner
  EXPECT_EQ(0,labelmove.getConnectedLabel());

}

TEST_F(TestMoveLabelAlongChain, checkAndApply)
{
  //IngredientsType ingredients;
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(false);
  ingredients.setPeriodicY(false);
  ingredients.setPeriodicZ(false);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(0,8,8);
  ingredients.modifyMolecules().addMonomer(2,8,8);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  ingredients.modifyMolecules().addMonomer(4,8,8);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  ingredients.modifyMolecules().addMonomer(6,8,8);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  ingredients.modifyMolecules().addMonomer(6,6,7);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  ingredients.modifyMolecules().addMonomer(4,6,7);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  ingredients.modifyMolecules().addMonomer(2,6,7);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  ingredients.modifyMolecules().addMonomer(0,6,7);
  ingredients.modifyMolecules().connect(ingredients.getMolecules().size()-2,ingredients.getMolecules().size()-1);
  
  ingredients.setNumTendomers(1);
  ingredients.setNumCrossLinkers(0);
  ingredients.setNumLabelsPerTendomerArm(1);
  ingredients.modifyMolecules()[3].setLabel(3);
  ingredients.modifyMolecules()[4].setLabel(4);
  ingredients.setNumMonomersPerChain(4);

  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(3,ingredients.modifyMolecules()[3].getLabel());
  EXPECT_EQ(4,ingredients.modifyMolecules()[4].getLabel());

  MoveLabelAlongChain labelmove;

  // init: not labeled monomers and thus should be rejected
  labelmove.init(ingredients,0);
  EXPECT_FALSE(labelmove.check(ingredients));
  labelmove.init(ingredients,1);
  EXPECT_FALSE(labelmove.check(ingredients));
  labelmove.init(ingredients,5);
  EXPECT_FALSE(labelmove.check(ingredients));
  labelmove.init(ingredients,6);
  EXPECT_FALSE(labelmove.check(ingredients));
  
  //labeled monomers hitting the 'wall'/chain ends
  labelmove.init(ingredients,3,MoveLabelAlongChain::RIGHT);
  EXPECT_FALSE(labelmove.check(ingredients));
  labelmove.init(ingredients,4,MoveLabelAlongChain::LEFT);
  EXPECT_FALSE(labelmove.check(ingredients));
  //labeled monomers move to non occupied place 
  labelmove.init(ingredients,3,MoveLabelAlongChain::LEFT);
  EXPECT_TRUE(labelmove.check(ingredients));
  labelmove.init(ingredients,4,MoveLabelAlongChain::RIGHT);
  EXPECT_TRUE(labelmove.check(ingredients));
  
  //
  labelmove.init(ingredients,4,MoveLabelAlongChain::RIGHT);
  labelmove.apply(ingredients);
  EXPECT_EQ(0,ingredients.getMolecules()[4].getLabel());
  EXPECT_EQ(4,ingredients.getMolecules()[5].getLabel());
  labelmove.init(ingredients,5,MoveLabelAlongChain::LEFT);
  EXPECT_TRUE(labelmove.check(ingredients));
  labelmove.apply(ingredients);
  //
  ingredients.modifyMolecules()[5].setLabel(4);
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(4,ingredients.modifyMolecules()[5].getLabel());
  //put another label and check if they see each other and thus
  //reject the move
  labelmove.init(ingredients,4,MoveLabelAlongChain::RIGHT);
  EXPECT_FALSE(labelmove.check(ingredients));
  //move the fifth monomer to the right to get space for the fourth
  labelmove.init(ingredients,5,MoveLabelAlongChain::RIGHT);
  EXPECT_TRUE(labelmove.check(ingredients));
  labelmove.apply(ingredients);
  
  EXPECT_EQ(4,ingredients.getMolecules()[6].getLabel());
  EXPECT_EQ(0,ingredients.getMolecules()[5].getLabel());
  //move label and move bond 
  labelmove.init(ingredients,4,MoveLabelAlongChain::RIGHT);
  EXPECT_TRUE(labelmove.check(ingredients));
  labelmove.apply(ingredients);
  EXPECT_TRUE ( ingredients.getMolecules().areConnected(3,5) );
  EXPECT_TRUE ( ingredients.getMolecules().areConnected(4,5) );
  EXPECT_FALSE( ingredients.getMolecules().areConnected(3,4) );
  labelmove.init(ingredients,3,MoveLabelAlongChain::LEFT);
  EXPECT_TRUE(labelmove.check(ingredients));
  labelmove.apply(ingredients);
  EXPECT_TRUE ( ingredients.getMolecules().areConnected(2,5) );
  EXPECT_TRUE ( ingredients.getMolecules().areConnected(2,3) );

}
