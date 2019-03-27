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
 * @brief Tests for UpdaterAddLinearChains
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
#include <LeMonADE/updater/UpdaterAddLinearChains.h>

using namespace std;




/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class TestUpdaterAddLinearChains: public ::testing::Test{
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

TEST_F(TestUpdaterAddLinearChains, Constructor)
{
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo <uint8_t> >,FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;

  IngredientsType ingredients;

  //Constructor call
  UpdaterAddLinearChains<IngredientsType> Timmy(ingredients,1,2);
  //check initialisation list
  EXPECT_FALSE(Timmy.getIsSolvent());
  EXPECT_EQ(2,Timmy.getNMonomerPerChain());
  EXPECT_EQ(1,Timmy.getNChain());
  EXPECT_DOUBLE_EQ(0.0,Timmy.getDensity());
}

TEST_F(TestUpdaterAddLinearChains, TestUpdater)
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

  UpdaterAddLinearChains<IngredientsType> Tommy(ingredients, 4, 16);

   // first execution
  EXPECT_NO_THROW(Tommy.initialize());
  //repeat execution, should not do anything
  EXPECT_TRUE(Tommy.execute());

  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ((4*16),ingredients.getMolecules().size());
  // check (default) monomer taggs
  EXPECT_EQ(1,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(1,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_EQ(2,ingredients.getMolecules()[15].getAttributeTag());

  // check connectivity
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(0));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(1));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(14));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(15));

  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(16));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(17));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(30));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(31));

}

TEST_F(TestUpdaterAddLinearChains, TestCompressedSolvent)
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

  // add some chains longer than 1 but allow solvent compression
  UpdaterAddLinearChains<IngredientsType> Tommy_1(ingredients, 4, 16, 3, 4, true);

  // first execution
  EXPECT_NO_THROW(Tommy_1.initialize());
  //repeat execution, should not do anything
  EXPECT_TRUE(Tommy_1.execute());
  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ((4*16),ingredients.getMolecules().size());
  // check explicitly set  monomer taggs
  EXPECT_EQ(3,ingredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(3,ingredients.getMolecules()[2].getAttributeTag());
  EXPECT_EQ(4,ingredients.getMolecules()[15].getAttributeTag());

  // expect no compressed solvent monomers (no chains of length 1)
  EXPECT_EQ(0,ingredients.getCompressedOutputIndices().size());

  // add some single monomers to same ingredients and try to compress
  UpdaterAddLinearChains<IngredientsType> Tommy_2(ingredients, 16, 1, 5, 5, true);
  EXPECT_NO_THROW(Tommy_2.initialize());
  EXPECT_TRUE(Tommy_2.execute());

  // use tags to check the new monomers
  EXPECT_EQ((5*16),ingredients.getMolecules().size());
  EXPECT_EQ(5,ingredients.getMolecules()[64].getAttributeTag());
  EXPECT_EQ(5,ingredients.getMolecules()[79].getAttributeTag());

  // get map of indizees size and index range
  EXPECT_EQ(1,ingredients.getCompressedOutputIndices().size());
  EXPECT_EQ(64,ingredients.getCompressedOutputIndices().begin()->first);
  EXPECT_EQ(79,ingredients.getCompressedOutputIndices().begin()->second);

  // add some more single monomers to same ingredients but do not compress
  UpdaterAddLinearChains<IngredientsType> Tommy_3(ingredients, 16, 1, 6, 6, false);
  EXPECT_NO_THROW(Tommy_3.initialize());
  EXPECT_TRUE(Tommy_3.execute());

  // use tags to check the new monomers
  EXPECT_EQ((6*16),ingredients.getMolecules().size());
  EXPECT_EQ(6,ingredients.getMolecules()[80].getAttributeTag());
  EXPECT_EQ(6,ingredients.getMolecules()[95].getAttributeTag());

  // same number of compressed monomers as before
  EXPECT_EQ(1,ingredients.getCompressedOutputIndices().size());
  EXPECT_EQ(64,ingredients.getCompressedOutputIndices().begin()->first);
  EXPECT_EQ(79,ingredients.getCompressedOutputIndices().begin()->second);
}

