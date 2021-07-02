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

#include <iostream>
#include <cstdio>

#include "gtest/gtest.h"

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureSystemInformationDendrimer.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>

class TestFeatureSystemInformationDendrimer: public ::testing::Test{
public:
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


TEST(TestFeatureSystemInformationDendrimer,ReadWrite)
{

    typedef LOKI_TYPELIST_1(FeatureSystemInformationDendrimer) Features;
    typedef ConfigureSystem<VectorInt3,Features> Config;
    typedef Ingredients<Config> IngredientsType;

    IngredientsType ingredients;

    constexpr uint32_t generation(5);
	  constexpr uint32_t spacerlength(2);
	  constexpr uint32_t coreFunctionality(4);
	  constexpr uint32_t branchingPointFunctionality(3);
    
    ingredients.setGeneration(generation);
    ingredients.setSpacerLength(spacerlength);
    ingredients.setCoreFunctionality(coreFunctionality);
    ingredients.setBranchingPointFunctionality(branchingPointFunctionality);

    EXPECT_NO_THROW(ingredients.synchronize());

    // Check getter and setter
    EXPECT_TRUE(ingredients.getGeneration() == 5);
    EXPECT_TRUE(ingredients.getSpacerLength() == 2);
    EXPECT_TRUE(ingredients.getCoreFunctionality() == 4);
    EXPECT_TRUE(ingredients.getBranchingPointFunctionality() == 3);

    // Check writing/reading to the bfm file
    AnalyzerWriteBfmFile<IngredientsType> bfmWriter("TestFeatureSystemInformationDendrimer.bfm", ingredients);
    bfmWriter.initialize();
    bfmWriter.execute();
    bfmWriter.cleanup();
    
    IngredientsType ingredientsRead;

    UpdaterReadBfmFile<IngredientsType> bfmReader("TestFeatureSystemInformationDendrimer.bfm",ingredientsRead,UpdaterReadBfmFile<IngredientsType>::READ_STEPWISE);
    bfmReader.initialize();

    EXPECT_TRUE(ingredients.getGeneration() == 5);
    EXPECT_TRUE(ingredients.getSpacerLength() == 2);
    EXPECT_TRUE(ingredients.getCoreFunctionality() == 4);
    EXPECT_TRUE(ingredients.getBranchingPointFunctionality() == 3);

    EXPECT_EQ(0,std::remove("TestFeatureSystemInformationDendrimer.bfm"));

}
