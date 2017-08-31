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
 * @brief Tests for the class FeatureFixedMonomers
 *
 * @author Martin
 * @date 07.07.2014
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <iostream>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

using namespace std;

// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class TestFeatureFixedMonomers
 * @brief checking different classes and functions depending to FeatureFixedMonomers
 * */
/*****************************************************************************/
class TestFeatureFixedMonomers: public ::testing::Test{
protected:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureFixedMonomers) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> MyIngredients;

public:

  //dummy move class used to check response to unknown move type
  class UnknownMove:public MoveBase
  {
	  public:
		template<class IngredientsType> bool check(IngredientsType& ingredients) const
		{
		return ingredients.checkMove(ingredients,*this);
		}

		template<class IngredientsType> void apply(IngredientsType& ingredients)
		{
		ingredients.applyMove(ingredients,*this);
		}

		template <class IngredientsType> void init(const IngredientsType& ingredients){};
  };

  /* suppress cout output for better readability -->un-/comment here:*/
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
  /* ** */
};

/*****************************************************************************/
TEST_F(TestFeatureFixedMonomers, CheckMoveableTag)
{
  MonomerMovableTag testTag;
  EXPECT_TRUE(testTag.getMovableTag());
  testTag.setMovableTag(false);
  EXPECT_FALSE(testTag.getMovableTag());

}

/*****************************************************************************/
TEST_F(TestFeatureFixedMonomers, CheckReadWrite)
{
  MyIngredients ingredients;

  ingredients.setBoxX(64);
  ingredients.setBoxY(64);
  ingredients.setBoxZ(64);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);

  ingredients.modifyMolecules().addMonomer(1,1,1);
  EXPECT_TRUE(ingredients.getMolecules()[0].getMovableTag());
  ingredients.modifyMolecules()[0].setMovableTag(false);
  EXPECT_FALSE(ingredients.getMolecules()[0].getMovableTag());

  ingredients.modifyMolecules().addMonomer(2,3,1);
  ingredients.modifyMolecules().addMonomer(3,5,1);
  ingredients.modifyMolecules().addMonomer(4,7,1);

  string filename("tests/tmpTestFeatureFixedMonomers.bfm");
  AnalyzerWriteBfmFile<MyIngredients> writeobject(filename, ingredients);

  writeobject.initialize();

  writeobject.execute();

  MyIngredients ingredients2;
  FileImport<MyIngredients> fileImport(filename,ingredients2);

  //scan file for !mcs and read-in first frame
  fileImport.initialize();


  EXPECT_FALSE(ingredients2.getMolecules()[0].getMovableTag());
  EXPECT_TRUE(ingredients2.getMolecules()[1].getMovableTag());
  EXPECT_TRUE(ingredients2.getMolecules()[2].getMovableTag());
  EXPECT_TRUE(ingredients2.getMolecules()[3].getMovableTag());
  remove(filename.c_str());

  MyIngredients ingredients3;
  ingredients3.modifyMolecules().addMonomer(1,1,1);
  ReadFixedMonomers <MyIngredients> command(ingredients3);

  //check reaction to correct standard input
  istringstream inputStandard("\n1-1:1");
  command.setInputStream(&inputStandard);
  EXPECT_NO_THROW(command.execute());
  EXPECT_FALSE(ingredients3.getMolecules()[0].getMovableTag());

  //check reaction to spacers
  istringstream inputSpacer("\n  1  - 1 :    1   ");
  command.setInputStream(&inputSpacer);
  EXPECT_NO_THROW(command.execute());
  EXPECT_FALSE(ingredients3.getMolecules()[0].getMovableTag());

  //check reaction to empty input
  istringstream inputEmpty("\n ");
  command.setInputStream(&inputEmpty);
  EXPECT_THROW(command.execute(),std::runtime_error);

  //check reaction to wrong format
  istringstream inputWrong1("\na-b:c");
  command.setInputStream(&inputWrong1);
  EXPECT_THROW(command.execute(),std::runtime_error);

  //check reaction to no/wrong seperators(-, :)
  istringstream inputWrong2("\n1:1-2");
  command.setInputStream(&inputWrong2);
  EXPECT_THROW(command.execute(),std::runtime_error);

  istringstream inputWrong3("\n1-1-2");
  command.setInputStream(&inputWrong3);
  EXPECT_THROW(command.execute(),std::runtime_error);
}
/*****************************************************************************/

TEST_F(TestFeatureFixedMonomers, CheckMove)
{
  MyIngredients ingredients;
  ingredients.modifyMolecules().addMonomer(1,1,1);

  UnknownMove move;
  EXPECT_TRUE(move.check(ingredients));

  MoveLocalSc localmove;
  localmove.init(ingredients);
  EXPECT_TRUE(localmove.check(ingredients));
  ingredients.modifyMolecules()[0].setMovableTag(false);
  EXPECT_FALSE(localmove.check(ingredients));
}
