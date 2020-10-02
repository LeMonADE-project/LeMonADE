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
 * @brief Tests for UpdaterAbstractCreate
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/UpdaterSwellBox.h>
#include <LeMonADE/feature/FeatureAttributes.h>

class TestUpdaterSwellBox: public ::testing::Test{
public:

  typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;
  
  IngredientsType ingredients;
  typename IngredientsType::molecules_type& setMolies = ingredients.modifyMolecules(); 
  const typename IngredientsType::molecules_type& getMolies = ingredients.getMolecules(); 
  //redirect cout output
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };
  void setConfig()
  {
    ingredients.setBoxX(16);
    ingredients.setBoxY(16);
    ingredients.setBoxZ(16);
    ingredients.setPeriodicX(false);
    ingredients.setPeriodicY(false);
    ingredients.setPeriodicZ(false);
    ingredients.modifyBondset().addBFMclassicBondset();
    EXPECT_NO_THROW(ingredients.synchronize());
    constexpr uint32_t nMonomers(10);
    ingredients.modifyMolecules().resize(nMonomers);
    uint32_t ID(0);
    VectorInt3 COM(8,8,8);
    VectorInt3 Bond(0,0,0);
    setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates( 2, 0, 0);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates(-2, 0, 0);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates( 0, 2, 0);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates( 0,-2, 0);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates( 0, 0, 2);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates( 0, 0,-2);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    for (auto i=1; i<7;i++) setMolies.connect(0,i);
    Bond.setAllCoordinates( 0, 0,-8);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    Bond.setAllCoordinates( 2, 0,-8);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
    setMolies.connect(7,8);
    Bond.setAllCoordinates( 6, 0,-8);setMolies[ID].modifyVector3D()=COM+Bond;ID++;
  }
private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
  

};

TEST_F(TestUpdaterSwellBox, CheckPeriodicityCheck)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(false);
  ingredients.setPeriodicZ(false);
  ingredients.modifyBondset().addBFMclassicBondset();
  EXPECT_NO_THROW(ingredients.synchronize());
  UpdaterSwellBox<IngredientsType> Guenther(ingredients,128, 32);
  EXPECT_ANY_THROW(Guenther.initialize());
}
TEST_F(TestUpdaterSwellBox, ClusterDetection)
{
  setConfig();
  EXPECT_NO_THROW(ingredients.synchronize());
  UpdaterSwellBox<IngredientsType> Guenther(ingredients,128, 32);
  EXPECT_NO_THROW(Guenther.initialize());
  EXPECT_EQ(3,Guenther.getNMolecules());
  EXPECT_EQ(7,Guenther.getLargestClusterSize());
}

TEST_F(TestUpdaterSwellBox, IncreaseBoxSizeClusterInRange)
{
  setConfig();
  EXPECT_NO_THROW(ingredients.synchronize());
  UpdaterSwellBox<IngredientsType> Heiner(ingredients,128, 32,7);
  EXPECT_NO_THROW(Heiner.initialize());
  EXPECT_EQ(48,ingredients.getBoxX());
  EXPECT_EQ(48,ingredients.getBoxY());
  EXPECT_EQ(48,ingredients.getBoxZ());
}
TEST_F(TestUpdaterSwellBox, IncreaseBoxSizeClusterOutOfRange)
{
  setConfig();
  EXPECT_NO_THROW(ingredients.synchronize());
  UpdaterSwellBox<IngredientsType> Heiner(ingredients,128, 32,6);
  EXPECT_NO_THROW(Heiner.initialize());
  EXPECT_EQ(16,ingredients.getBoxX());
  EXPECT_EQ(16,ingredients.getBoxY());
  EXPECT_EQ(16,ingredients.getBoxZ());
}
