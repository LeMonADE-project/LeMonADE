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
 * @brief Tests for the class FeatureJumps
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
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/FeatureJumps.h>

using namespace std;
// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class TestFeatureJumps
 * @brief checking different classes and functions depending to FeatureFixedMonomers
 * */
/*****************************************************************************/
class TestFeatureJumps: public ::testing::Test{
protected:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureJumps) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> MyIngredients;
  MyIngredients ingredients;

  /* suppress cout output for better readability -->un-/comment here:*/
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
  /* ** */
};


/*****************************************************************************/
TEST_F(TestFeatureJumps, CheckWrite)
{
  //set box size
  ingredients.setBoxX(32);
  ingredients.setBoxY(32);
  ingredients.setBoxZ(32);
  //set box periodicity
  ingredients.setPeriodicX(false);
  ingredients.setPeriodicY(false);
  ingredients.setPeriodicZ(true);

  //add bond to the bondset
  //add mirrored bond vector for consistency
  ingredients.modifyBondset().addBond(1,1,2,48);
  ingredients.modifyBondset().addBond(-1,-1,-2,49); 
  uint i = 0;
  for(i=0;i<16;i++){
    ingredients.modifyMolecules().addMonomer(VectorInt3(i,i,2*i+32));
  }
  for(i=0;i<15;i++){
    ingredients.modifyMolecules().connect(i,i+1);
  }

  string filename("tests/tmpTestFeatureJumps.bfm");
  AnalyzerWriteBfmFile<MyIngredients> writeobject(filename, ingredients);
  writeobject.initialize();
  writeobject.execute();
  EXPECT_EQ(0,remove("tests/tmpTestFeatureJumps.bfm"));
}


/*****************************************************************************/
TEST_F(TestFeatureJumps, CheckRead)
{


}
