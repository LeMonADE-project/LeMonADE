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
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/utility/Vector3D.h>

#include <LeMonADE/updater/moves/MoveBreak.h>

class TestFeatureConnection : public ::testing::Test
{
public:
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureConnectionSc,FeatureExcludedVolumeSc<> ) Features;
  typedef ConfigureSystem<VectorInt3,Features,3> Config;
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
void setConfig()
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
    ingredients.synchronize(ingredients);
}
private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
};

TEST_F(TestFeatureConnection,LatticeOccupation)
{
  setConfig();
  EXPECT_EQ(0,ingredients.getIdFromLattice(0,0,0));    
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(),ingredients.getIdFromLattice(1,0,0));  
  EXPECT_EQ(std::numeric_limits<uint32_t>::max(),ingredients.getIdFromLattice(1,17,32));  
  EXPECT_EQ(1,ingredients.getIdFromLattice(2,0,0));  
}

TEST_F(TestFeatureConnection,MoveAddMonomerSc)
{
   //prepare ingredients
    setConfig();
    //add monomer with move 
    MoveAddMonomerSc<> move;
    move.init(ingredients);
    VectorInt3 pos(0,2,0);
    move.setPosition(pos);
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(2,ingredients.getIdFromLattice(pos));
}

TEST_F(TestFeatureConnection,MoveLocalSc)
{
  //prepare ingredients
    setConfig();
    //move monomers and check the lattice occupation 
    MoveLocalSc move;

    VectorInt3 dir(0,1,0);
    move.init(ingredients,0,dir);
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(0,ingredients.getIdFromLattice(dir));
    EXPECT_EQ(std::numeric_limits<uint32_t>::max(),ingredients.getIdFromLattice(0,0,0));
    
    move.init(ingredients,1,dir);
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(1,ingredients.getIdFromLattice(2,1,0));
    EXPECT_EQ(std::numeric_limits<uint32_t>::max(),ingredients.getIdFromLattice(2,0,0));
}
TEST_F(TestFeatureConnection, MoveConnectSc)
{
  ingredients.setBoxX(16);
  ingredients.setBoxY(16);
  ingredients.setBoxZ(16);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.modifyMolecules().addMonomer(6,8,8);
  ingredients.modifyMolecules().addMonomer(10,8,8);
  ingredients.modifyMolecules().addMonomer(8,8,8);
  ingredients.modifyMolecules().addMonomer(8,6,8);
  ingredients.modifyMolecules().addMonomer(8,10,8);

  EXPECT_NO_THROW(ingredients.synchronize());
  EXPECT_EQ(5,ingredients.getMolecules().size());

  MoveConnectSc move;
  VectorInt3 dir(2,0,0);
  move.init(ingredients,2,dir);
  EXPECT_TRUE( move.check(ingredients));
  EXPECT_EQ(1,move.getPartner());
  move.apply(ingredients); //first connection 
  
  dir.setAllCoordinates(-2,0,0);
  move.init(ingredients,2,dir);
  EXPECT_TRUE( move.check(ingredients)); 
  move.apply(ingredients); //second connection 
  
  dir.setAllCoordinates(0,-2,0);
  move.init(ingredients,2,dir);
  EXPECT_TRUE(move.check(ingredients));
  move.apply(ingredients); //third connection 
  
  dir.setAllCoordinates(0,2,0);
  move.init(ingredients,2,dir );
  EXPECT_FALSE(move.check(ingredients));  //would be the fourth connection 
  
  dir.setAllCoordinates(0,-2,0);
  move.init(ingredients,2,dir);
  EXPECT_FALSE(move.check(ingredients)); //already connected 
 
}
