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


#include "gtest/gtest.h"
#include <cstdio>
#include <sstream>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBreak.h>
#include <LeMonADE/utility/Vector3D.h>

#include <LeMonADE/updater/moves/MoveBreak.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>

class TestFeatureBreak : public ::testing::Test
{
public:
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureBreak, FeatureExcludedVolumeSc<> ) Features;
  typedef ConfigureSystem<VectorInt3,Features,17> Config; 
  typedef Ingredients<Config> Ing;
  Ing ingredients;

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


TEST_F(TestFeatureBreak, MoveBreakBad)
{
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyMolecules().resize(2);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules().connect(0,1);
    ingredients.modifyBondset().addBond(2,0,0,77);
    ingredients.modifyBondset().addBond(-2,0,0,75);
    ingredients.synchronize(ingredients);
    
    MoveBreak breakbad;
    breakbad.init(ingredients);
    EXPECT_TRUE(breakbad.check(ingredients)); 
    breakbad.apply(ingredients);
    
    breakbad.init(ingredients,0,1);
    EXPECT_FALSE(breakbad.check(ingredients));
}
