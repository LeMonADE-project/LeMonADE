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
 * @brief Tests for the class FeatureBoltzmann
 *
 * @author Martin
 * @date 07.07.2014
 * */
/*****************************************************************************/

#include <iostream>

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>

using namespace std;
// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class FeatureBoltzmannTest
 * @brief checking probability based move acceptance
 * */
/*****************************************************************************/
class FeatureBoltzmannTest: public ::testing::Test{

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

  /* suppress cout output for better readability -->un/comment here:*/
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
/**
 * @fn TEST_F(FeatureBoltzmannTest, CheckMove)
 * @brief Test function for CheckMove function of feature Boltzmann
 * */
/*****************************************************************************/
TEST_F(FeatureBoltzmannTest, CheckMove)
{
  UnknownMove testmove;
  FeatureBoltzmann feature;
  //setting move probability to one and try several times if checkMove returns true as expected
  testmove.resetProbability();
  EXPECT_TRUE(testmove.check(feature));
  EXPECT_TRUE(testmove.check(feature));
  EXPECT_TRUE(testmove.check(feature));

  //setting move probability to zero and try several times if checkMove returns false as expected
  testmove.multiplyProbability(0.0);
  EXPECT_FALSE(testmove.check(feature));
  EXPECT_FALSE(testmove.check(feature));
  EXPECT_FALSE(testmove.check(feature));
}
