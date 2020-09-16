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
#include <LeMonADE/feature/FeatureLabel.h>
#include <LeMonADE/utility/Vector3D.h>

class TestFeatureLabel : public ::testing::Test
{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureLabel ) Features;
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


TEST_F(TestFeatureLabel,MonomerLabel)
{
  //prepare ingredients
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyMolecules().resize(1);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    EXPECT_EQ (0, ingredients.getMolecules()[0].getLabel() );
    ingredients.modifyMolecules()[0].setLabel(3);
    EXPECT_EQ (3, ingredients.getMolecules()[0].getLabel() );
    ingredients.modifyMolecules()[0].setLabel(17);
    EXPECT_EQ (17, ingredients.getMolecules()[0].getLabel() );
}
TEST_F(TestFeatureLabel,LabelPartner)
{
  //prepare ingredients
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.modifyMolecules().resize(3);
    ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
    ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
    ingredients.modifyMolecules()[2].setAllCoordinates(4,0,0);
    ingredients.modifyMolecules()[0].setLabel(3);
    ingredients.modifyMolecules()[1].setLabel(4);
    ingredients.modifyMolecules().connect(0,1);
    ingredients.modifyMolecules().connect(1,2);
    ingredients.synchronize();
    EXPECT_EQ(2,ingredients.getLabelPartner(0));
    EXPECT_EQ(1,ingredients.getLabelPartner(1));
    EXPECT_EQ(0,ingredients.getLabelPartner(2));
    ingredients.modifyMolecules()[1].setLabel(3);
    ingredients.synchronize();
    EXPECT_EQ(0,ingredients.getLabelPartner(0));
    EXPECT_EQ(0,ingredients.getLabelPartner(1));
    EXPECT_EQ(0,ingredients.getLabelPartner(2));
    ingredients.modifyMolecules()[0].setLabel(0);
    ingredients.modifyMolecules()[1].setLabel(2);
    ingredients.modifyMolecules()[2].setLabel(3);
    ingredients.synchronize();
    EXPECT_EQ(3,ingredients.getLabelPartner(1));
    EXPECT_EQ(2,ingredients.getLabelPartner(2));
    
}


TEST_F(TestFeatureLabel,ReadLabelPartner)
{
  //this test is for testing the reaction of the read class to incorrectly
  //formatted input
  Ing ingredients;
  ReadLabel<Ing> read(ingredients);

  ingredients.modifyMolecules().resize(10);

  std::stringstream stream1;
  read.setInputStream(&stream1);
  stream1<<"\ni am a parrot\n";
  EXPECT_THROW(read.execute(),std::runtime_error);


  std::stringstream stream2;
  read.setInputStream(&stream2);
  stream2<<"\n00011:2:11110:3\n";
  EXPECT_THROW(read.execute(),std::runtime_error);


  std::stringstream stream3;
  read.setInputStream(&stream3);
  stream3<<"\n00011:2:11110-3\n";
  EXPECT_THROW(read.execute(),std::runtime_error);


  std::stringstream stream4;
  read.setInputStream(&stream4);
  stream4<<"\n00011:2:11230-3\n";
  EXPECT_THROW(read.execute(),std::runtime_error);

  std::stringstream stream6;
  read.setInputStream(&stream6);
  stream6<<"\n00011:2-11230:3\n";
  EXPECT_THROW(read.execute(),std::runtime_error);
  
  std::stringstream stream5;
  read.setInputStream(&stream5);
  stream5<<"\n00011:2-11110:3\n";
  EXPECT_NO_THROW(read.execute());

  for(size_t i = 0; i < 3; i++)
    EXPECT_EQ(0,ingredients.getMolecules()[i].getLabel());

  for(size_t i = 3; i < 9; i++)
    EXPECT_TRUE(ingredients.getMolecules()[i].getLabel()>0);

  for(size_t i = 3; i < 5; i++)
    EXPECT_EQ(2,ingredients.getMolecules()[i].getLabel());

  for(size_t i = 5; i < 9; i++)
    EXPECT_EQ(3,ingredients.getMolecules()[i].getLabel());

  for(size_t i = 9; i < 10; i++)
    EXPECT_EQ(0,ingredients.getMolecules()[i].getLabel());

}
TEST_F(TestFeatureLabel, WriteLabelPartner)
{
  Ing ingredients;
  WriteLabel<Ing> write(ingredients);
  ingredients.modifyMolecules().resize(10);
  std::stringstream writeStream;
  ingredients.setBoxX(12);
  ingredients.setBoxY(12);
  ingredients.setBoxZ(12);
  ingredients.modifyMolecules().setAge(10000);
  ingredients.setPeriodicX(1);
  ingredients.setPeriodicY(1);
  ingredients.setPeriodicZ(1);
  ingredients.setNumMonomersPerChain(5);
  ingredients.setNumTendomers(1);
  ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
  ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
  ingredients.modifyMolecules()[2].setAllCoordinates(2,2,0);
  ingredients.modifyMolecules()[3].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[4].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[5].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[6].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[7].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[8].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[9].setAllCoordinates(2,0,2);
  ingredients.modifyMolecules()[0].setLabel(2);
  ingredients.modifyMolecules()[1].setLabel(2);
  ingredients.modifyMolecules()[2].setLabel(2);
  ingredients.modifyMolecules()[3].setLabel(0);
  ingredients.modifyMolecules()[4].setLabel(0);
  ingredients.modifyMolecules()[5].setLabel(4);
  ingredients.modifyMolecules()[6].setLabel(0);
  ingredients.modifyMolecules()[7].setLabel(4);
  ingredients.modifyMolecules()[8].setLabel(0);
  ingredients.modifyMolecules()[9].setLabel(4);
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
  EXPECT_EQ(checkIngredients.getMolecules().size(),10);
  EXPECT_EQ(2,checkIngredients.getMolecules()[0].getLabel());
  EXPECT_EQ(2,checkIngredients.getMolecules()[1].getLabel());
  EXPECT_EQ(2,checkIngredients.getMolecules()[2].getLabel());
  EXPECT_EQ(0,checkIngredients.getMolecules()[3].getLabel());
  EXPECT_EQ(0,checkIngredients.getMolecules()[4].getLabel());
  EXPECT_EQ(4,checkIngredients.getMolecules()[5].getLabel());
  EXPECT_EQ(0,checkIngredients.getMolecules()[6].getLabel());
  EXPECT_EQ(4,checkIngredients.getMolecules()[7].getLabel());
  EXPECT_EQ(0,checkIngredients.getMolecules()[8].getLabel());
  EXPECT_EQ(4,checkIngredients.getMolecules()[9].getLabel());
  //remove temporary file tmp.bfm (and let test fail if not removed)
  EXPECT_EQ(0,remove(filename.c_str()));    
}

TEST_F(TestFeatureLabel,MoveAddMonomerSc)
{
  //prepare ingredients
    ingredients.setBoxX(12);
    ingredients.setBoxY(12);
    ingredients.setBoxZ(12);
    ingredients.setPeriodicX(1);
    ingredients.setPeriodicY(1);
    ingredients.setPeriodicZ(1);
    ingredients.modifyBondset().addBFMclassicBondset();
    ingredients.synchronize(ingredients);
    
    ingredients.setNumTendomers(1);
    ingredients.setNumCrossLinkers(0);
    ingredients.setNumLabelsPerTendomerArm(1);

    ingredients.setNumMonomersPerChain(3);
    ingredients.synchronize(ingredients);
    MoveAddMonomerSc<> move;
    MoveConnectSc connection;
    move.init(ingredients);
    VectorInt3 pos(0,2,0);
    move.setPosition(pos);
    move.setLabel(0);
    EXPECT_EQ(0,move.getLabel());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    
    move.init(ingredients);
    pos.setAllCoordinates(0,4,0);
    move.setPosition(pos);
    move.setLabel(0);
    EXPECT_EQ(0,move.getLabel());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    connection.init(ingredients,ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
    connection.check(ingredients);
    connection.apply(ingredients);
    
    move.init(ingredients);
    pos.setAllCoordinates(0,6,0);
    move.setPosition(pos);
    move.setLabel(4);
    EXPECT_EQ(4,move.getLabel());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(4,ingredients.getMolecules()[2].getLabel());
    connection.init(ingredients,ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
    connection.check(ingredients);
    connection.apply(ingredients);
    
    move.init(ingredients);
    pos.setAllCoordinates(2,6,1);
    move.setPosition(pos);
    move.setLabel(5);
    EXPECT_EQ(5,move.getLabel());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    EXPECT_EQ(5,ingredients.getMolecules()[3].getLabel());
    connection.init(ingredients,ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
    connection.check(ingredients);
    connection.apply(ingredients);

    
    move.init(ingredients);
    pos.setAllCoordinates(2,4,1);
    move.setPosition(pos);
    move.setLabel(0);
    EXPECT_EQ(0,move.getLabel());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    connection.init(ingredients,ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
    connection.check(ingredients);
    connection.apply(ingredients);

    
    move.init(ingredients);
    pos.setAllCoordinates(2,2,1);
    move.setPosition(pos);
    move.setLabel(0);
    EXPECT_EQ(0,move.getLabel());
    EXPECT_TRUE(move.check(ingredients));
    move.apply(ingredients);
    connection.init(ingredients,ingredients.getMolecules().size()-2, ingredients.getMolecules().size()-1);
    connection.check(ingredients);
    connection.apply(ingredients);

//     ingredients.synchronize();
    MoveLabelAlongChain lbmove;
    lbmove.init(ingredients,2,MoveLabelAlongChain::RIGHT);
    EXPECT_EQ(4, ingredients.getMolecules()[2].getLabel());
    EXPECT_FALSE(lbmove.check(ingredients));
    EXPECT_EQ(4,lbmove.getConnectedLabel());
    
    lbmove.init(ingredients,2,MoveLabelAlongChain::LEFT);
    EXPECT_EQ(4, ingredients.getMolecules()[2].getLabel());
    EXPECT_TRUE(lbmove.check(ingredients));
    EXPECT_EQ(4,lbmove.getConnectedLabel());
    
}

