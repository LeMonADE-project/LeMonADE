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

class TestFeatureConnectionSc : public ::testing::Test
{
public:
  typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureConnectionSc,FeatureExcludedVolumeSc<> ) Features;
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
    ingredients.modifyMolecules()[0].setNumMaxLinks(1);
    ingredients.modifyMolecules()[1].setNumMaxLinks(17);
    ingredients.modifyMolecules()[2].setNumMaxLinks(3);
    EXPECT_TRUE (ingredients.getMolecules()[0].isReactive() );
    EXPECT_TRUE (ingredients.getMolecules()[1].isReactive() );
    EXPECT_FALSE(ingredients.getMolecules()[2].isReactive() );
    EXPECT_EQ(ingredients.getMolecules()[0].getNumMaxLinks(),1);
    EXPECT_EQ(ingredients.getMolecules()[1].getNumMaxLinks(),17);
    EXPECT_EQ(ingredients.getMolecules()[2].getNumMaxLinks(),3);
    ingredients.synchronize(ingredients);
    
    size_t idx = ingredients.modifyMolecules().addMonomer(34,7,8);
    ingredients.modifyMolecules()[idx].setReactive(true);
    ingredients.modifyMolecules()[idx].setNumMaxLinks(1);

    ingredients.synchronize(ingredients);

    EXPECT_TRUE  (ingredients.getMolecules()[0].isReactive() == ingredients.getMolecules()[0].isReactive());
    EXPECT_FALSE (ingredients.getMolecules()[2].isReactive() == ingredients.getMolecules()[0].isReactive());

    MonomerReactivity reactivity0=ingredients.getMolecules()[0].getMonomerReactivity();
    MonomerReactivity reactivity1=ingredients.getMolecules()[1].getMonomerReactivity();
    MonomerReactivity reactivity2=ingredients.getMolecules()[2].getMonomerReactivity();
    MonomerReactivity reactivity3=ingredients.getMolecules()[idx].getMonomerReactivity();

    EXPECT_TRUE  (reactivity0 == reactivity3);
    EXPECT_FALSE (reactivity0 == reactivity1);
    EXPECT_FALSE (reactivity0 == reactivity2);

    EXPECT_FALSE  (reactivity0 != reactivity3);
    EXPECT_TRUE (reactivity0 != reactivity1);
    EXPECT_TRUE (reactivity0 != reactivity2);
}

TEST_F(TestFeatureConnectionSc,ReadReactivity)
{
  //this test is for testing the reaction of the read class to incorrectly
  //formatted input
  Ing ingredients;
  ReadReactivity<Ing> read(ingredients);

  ingredients.modifyMolecules().resize(10);

  std::stringstream stream1;
  read.setInputStream(&stream1);
  stream1<<"\ni am a parrot\n";
  EXPECT_THROW(read.execute(),std::runtime_error);


  std::stringstream stream2;
  read.setInputStream(&stream2);
  stream2<<"\n1:2:2\n";
  EXPECT_THROW(read.execute(),std::runtime_error);


  std::stringstream stream3;
  read.setInputStream(&stream3);
  stream3<<"\n1-2-3\n";
  EXPECT_THROW(read.execute(),std::runtime_error);


  std::stringstream stream4;
  read.setInputStream(&stream4);
  stream4<<"\n1-a:4\n";
  EXPECT_THROW(read.execute(),std::runtime_error);

  std::stringstream stream5;
  read.setInputStream(&stream5);
  stream5<<"\n1-2:0/6\n";
  stream5<<"3-5:1/3\n";
  stream5<<"6-10:0/9\n";
  EXPECT_NO_THROW(read.execute());

  for(size_t i = 0; i < 2; i++)
      EXPECT_FALSE(ingredients.getMolecules()[i].isReactive());

  for(size_t i = 2; i < 5; i++)
	  EXPECT_TRUE(ingredients.getMolecules()[i].isReactive());

  for(size_t i = 5; i < ingredients.getMolecules().size(); i++)
  	  EXPECT_FALSE(ingredients.getMolecules()[i].isReactive());

  for(size_t i = 0; i < 2; i++)
	  EXPECT_EQ(ingredients.getMolecules()[i].getNumMaxLinks(),6);

  for(size_t i = 2; i < 5; i++)
	  EXPECT_EQ(ingredients.getMolecules()[i].getNumMaxLinks(),3);

  for(size_t i = 5; i < ingredients.getMolecules().size(); i++)
	  EXPECT_EQ(ingredients.getMolecules()[i].getNumMaxLinks(),9);
}
TEST_F(TestFeatureConnectionSc, WriteReactivity)
{
  Ing ingredients;
  WriteReactivity<Ing> write(ingredients);
  ingredients.modifyMolecules().resize(4);
  std::stringstream writeStream;
      ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.modifyMolecules().setAge(10000);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(2,2,0);
    ingredients.modifyMolecules()[3].setAllCoordinates(2,0,2);
    ingredients.modifyMolecules()[0].setReactive(true);
    ingredients.modifyMolecules()[1].setReactive(false);
    ingredients.modifyMolecules()[2].setReactive(false);
    ingredients.modifyMolecules()[3].setReactive(true);
    ingredients.modifyMolecules()[0].setNumMaxLinks(1);
    ingredients.modifyMolecules()[1].setNumMaxLinks(17);
    ingredients.modifyMolecules()[2].setNumMaxLinks(17);
    ingredients.modifyMolecules()[3].setNumMaxLinks(7);
    ingredients.synchronize(ingredients);
  
    //now create file with the information entered above
    std::string filename("tmpWriteReactivity.bfm");
    AnalyzerWriteBfmFile<Ing> outputFile(filename,ingredients,AnalyzerWriteBfmFile<Ing>::NEWFILE);
    outputFile.initialize();
    outputFile.execute();
    outputFile.closeFile();

    //now read the file back in and compare
    Ing checkIngredients;
    FileImport<Ing> inputFile(filename,checkIngredients);
    inputFile.initialize();

    //check age
    EXPECT_EQ(checkIngredients.getMolecules().getAge(),10000);
    //check number of monomers
    EXPECT_EQ(checkIngredients.getMolecules().size(),4);
    EXPECT_TRUE(checkIngredients.getMolecules()[0].isReactive()); 
    EXPECT_FALSE(checkIngredients.getMolecules()[1].isReactive()); 
    EXPECT_FALSE(checkIngredients.getMolecules()[2].isReactive()); 
    EXPECT_TRUE(checkIngredients.getMolecules()[3].isReactive()); 
    
    EXPECT_EQ(1,ingredients.getMolecules()[0].getNumMaxLinks());
    EXPECT_EQ(17,ingredients.getMolecules()[1].getNumMaxLinks());
    EXPECT_EQ(17,ingredients.getMolecules()[2].getNumMaxLinks());
    EXPECT_EQ(7,ingredients.getMolecules()[3].getNumMaxLinks());
    //remove temporary file tmp.bfm (and let test fail if not removed)
    EXPECT_EQ(0,remove(filename.c_str()));
     
}

TEST_F(TestFeatureConnectionSc,MoveAddMonomerSc)
{
  //prepare ingredients
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyMolecules().resize(2);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[0].setReactive(true);
    ingredients.modifyMolecules()[1].setReactive(false);
    ingredients.modifyMolecules()[0].setNumMaxLinks(1);
    ingredients.modifyMolecules()[1].setNumMaxLinks(17);
    EXPECT_TRUE (ingredients.getMolecules()[0].isReactive() );
    EXPECT_FALSE(ingredients.getMolecules()[1].isReactive() );
    EXPECT_EQ(ingredients.getMolecules()[0].getNumMaxLinks(),1);
    EXPECT_EQ(ingredients.getMolecules()[1].getNumMaxLinks(),17);
    ingredients.synchronize(ingredients);
    EXPECT_EQ(0,ingredients.getIdFromLattice(0,0,0));    
    MoveAddMonomerSc<> move;
    move.init(ingredients);
    VectorInt3 pos(0,2,0);
    move.setPosition(pos);
    move.setReactive(true);
    EXPECT_TRUE(move.isReactive());
    move.setNumMaxLinks(4);
    EXPECT_EQ(4,move.getNumMaxLinks());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(2,ingredients.getIdFromLattice(pos));
}

TEST_F(TestFeatureConnectionSc,MoveLocalSc)
{
  //prepare ingredients
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyMolecules().resize(2);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[0].setReactive(true);
    ingredients.modifyMolecules()[1].setReactive(false);
    ingredients.modifyMolecules()[0].setNumMaxLinks(1);
    ingredients.modifyMolecules()[1].setNumMaxLinks(17);
    EXPECT_TRUE (ingredients.getMolecules()[0].isReactive() );
    EXPECT_FALSE(ingredients.getMolecules()[1].isReactive() );
    EXPECT_EQ(ingredients.getMolecules()[0].getNumMaxLinks(),1);
    EXPECT_EQ(ingredients.getMolecules()[1].getNumMaxLinks(),17);
    ingredients.synchronize(ingredients);
    EXPECT_EQ(0,ingredients.getIdFromLattice(0,0,0));    
    EXPECT_EQ(std::numeric_limits<uint32_t>::max(),ingredients.getIdFromLattice(2,0,0));    
    MoveLocalSc move;

    VectorInt3 dir(0,1,0);
    move.init(ingredients,0,dir);
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(0,ingredients.getIdFromLattice(dir));
    
    move.init(ingredients,1,dir);
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(std::numeric_limits<uint32_t>::max(),ingredients.getIdFromLattice(2,1,0));
}
