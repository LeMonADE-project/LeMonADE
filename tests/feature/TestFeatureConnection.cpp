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
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/utility/Vector3D.h>

class TestFeatureConnectionSc : public ::testing::Test
{
public:
  typedef LOKI_TYPELIST_3(FeatureConnectionSc, FeatureBondset< >,FeatureBox ) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
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


TEST_F(TestFeatureConnectionSc,MonomerReactivitySetting)
{
  //prepare ingredients
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(1,0,30);
    ingredients.modifyMolecules()[0].setReactive(true);
    ingredients.modifyMolecules()[1].setReactive(true);
    ingredients.modifyMolecules()[2].setReactive(false);
    ingredients.modifyMolecules()[0].setNMaxBonds(1);
    ingredients.modifyMolecules()[1].setNMaxBonds(17);
    ingredients.modifyMolecules()[2].setNMaxBonds(3);
    EXPECT_TRUE (ingredients.getMolecules()[0].isReactive() );
    EXPECT_TRUE (ingredients.getMolecules()[1].isReactive() );
    EXPECT_FALSE(ingredients.getMolecules()[2].isReactive() );
    EXPECT_EQ(ingredients.getMolecules()[0].getNMaxBonds(),1);
    EXPECT_EQ(ingredients.getMolecules()[1].getNMaxBonds(),17);
    EXPECT_EQ(ingredients.getMolecules()[2].getNMaxBonds(),3);
    ingredients.synchronize(ingredients);
    
    size_t idx = ingredients.modifyMolecules().addMonomer(34,7,8);
    ingredients.modifyMolecules()[idx].setReactive(true);
    ingredients.modifyMolecules()[idx].setNMaxBonds(1);

    EXPECT_TRUE  (ingredients.getMolecules()[0].isReactive() == ingredients.getMolecules()[0].isReactive());
    EXPECT_FALSE (ingredients.getMolecules()[2].isReactive() == ingredients.getMolecules()[0].isReactive());

    MonomerReactivity reactivity0=ingredients.getMolecules()[0].getReactivity();
    MonomerReactivity reactivity1=ingredients.getMolecules()[1].getReactivity();
    MonomerReactivity reactivity2=ingredients.getMolecules()[2].getReactivity();
    MonomerReactivity reactivity3=ingredients.getMolecules()[idx].getReactivity();

    EXPECT_TRUE  (reactivity0 == reactivity3);
    EXPECT_FALSE (reactivity0 == reactivity1);
    EXPECT_FALSE (reactivity0 == reactivity2);

    EXPECT_FALSE  (reactivity0 != reactivity3);
    EXPECT_TRUE (reactivity0 != reactivity1);
    EXPECT_TRUE (reactivity0 != reactivity2);
}
// Test_F(TestFeatureConnectionSc,LatticeSetup)
// {
//   //prepare ingredients
//     ingredients.setBoxX(32);
//     ingredients.setBoxY(32);
//     ingredients.setBoxZ(32);
//     ingredients.setPeriodicX(1);
//     ingredients.setPeriodicY(1);
//     ingredients.setPeriodicZ(1);
//     ingredients.modifyMolecules().resize(3);
//     ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
//     ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
//     ingredients.modifyMolecules()[2].setAllCoordinates(1,0,30);
// 
//     ingredients.synchronize(ingredients);
//     
//     
// }
