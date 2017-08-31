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
 * @brief Tests for MonomerAttributeTag,FeatureAttribute
 * 
 * @date 30.06.2014
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up 
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class FeatureAttributesTest: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureAttributes) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> MyIngredients;
  
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
};

TEST_F(FeatureAttributesTest,MonomerAttributeTag)
{
  MonomerAttributeTag tag;
  EXPECT_EQ(0,tag.getAttributeTag());
  tag.setAttributeTag(-1);
  EXPECT_EQ(-1,tag.getAttributeTag());
  tag.setAttributeTag(15);
  EXPECT_EQ(15,tag.getAttributeTag());
}

TEST_F(FeatureAttributesTest,exportRead)
{
  MyIngredients ingredients;
  MyIngredients::molecules_type const& molecules= ingredients.getMolecules(); 
  FileImport<MyIngredients> file ("tests/attributesTest.test",ingredients);
  
  //scan file for !mcs and read-in first frame
  file.initialize();

  EXPECT_EQ(molecules[0].getAttributeTag(),0); /*this one has the default value*/
  EXPECT_EQ(molecules[1].getAttributeTag(),5);
  EXPECT_EQ(molecules[2].getAttributeTag(),5);
  EXPECT_EQ(molecules[3].getAttributeTag(),5);
  EXPECT_EQ(molecules[4].getAttributeTag(),2);
  EXPECT_EQ(molecules[6].getAttributeTag(),5);
  EXPECT_EQ(molecules[14].getAttributeTag(),5);
  EXPECT_EQ(molecules[15].getAttributeTag(),2);
  
  
}

TEST_F(FeatureAttributesTest,exportWrite)
{
  MyIngredients setupIngrediens;
  MyIngredients ingredients;
  
  MyIngredients::molecules_type const& setupMolecules= setupIngrediens.getMolecules(); 
  MyIngredients::molecules_type const& molecules= ingredients.getMolecules(); 
  
  setupIngrediens.setBoxX(100);
  setupIngrediens.setPeriodicX(true);
  setupIngrediens.setBoxY(100);
  setupIngrediens.setPeriodicY(true);
  setupIngrediens.setBoxZ(100);
  setupIngrediens.setPeriodicZ(true);
  setupIngrediens.modifyMolecules().resize(5);
  setupIngrediens.modifyMolecules()[0].setAttributeTag(1);
  setupIngrediens.modifyMolecules()[1].setAttributeTag(1);
  setupIngrediens.modifyMolecules()[4].setAttributeTag(2);
  //write to file and read back in
  AnalyzerWriteBfmFile<MyIngredients> outfile("tests/attributesTestOut.test",setupIngrediens,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
  outfile.initialize();

  FileImport<MyIngredients> infile ("tests/attributesTestOut.test",ingredients);
  
  //scan file for !mcs and read-in first frame
  infile.initialize();

  EXPECT_EQ(molecules[0].getAttributeTag(),1);
  EXPECT_EQ(molecules[1].getAttributeTag(),1);
  EXPECT_EQ(molecules[2].getAttributeTag(),0); /*has default value*/
  EXPECT_EQ(molecules[3].getAttributeTag(),0); /*has default value*/
  EXPECT_EQ(molecules[4].getAttributeTag(),2);
  
    //remove the temporary file
  EXPECT_EQ(0,remove("tests/attributesTestOut.test"));
}

TEST_F(FeatureAttributesTest,ReadAttributesClass)
{
  //this test is for testing the reaction of the read class to incorrectly 
  //formatted input
  MyIngredients ingredients;
  ReadAttributes<MyIngredients> read(ingredients);
  
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
  stream5<<"\n1-5:c\n";
  EXPECT_THROW(read.execute(),std::runtime_error);
  
  
}

TEST_F(FeatureAttributesTest,CopyConstructor)
{
	typedef MyIngredients::molecules_type MyMolecules;
	MyMolecules molecules1;
	
	//check if attributes are copied corriectly by copy constructor
	molecules1.resize(5);
	molecules1[0].setAttributeTag(1);
	molecules1[1].setAttributeTag(2);
	molecules1[2].setAttributeTag(3);
	molecules1[3].setAttributeTag(4);
	molecules1[4].setAttributeTag(5);
	//create new objects from molecules1
	
	typedef ConfigureSystem<VectorInt3,Features,8> Config8;
	typedef Ingredients<Config8> MyIngredients8;
	typedef MyIngredients8::molecules_type MyMolecules8;
	MyMolecules molecules2(molecules1);
	MyMolecules8 molecules3(molecules2);
	
	EXPECT_EQ(molecules1[0].getAttributeTag(),molecules2[0].getAttributeTag());
	EXPECT_EQ(molecules1[0].getAttributeTag(),molecules3[0].getAttributeTag());
	
	EXPECT_EQ(molecules1[1].getAttributeTag(),molecules2[1].getAttributeTag());
	EXPECT_EQ(molecules1[1].getAttributeTag(),molecules3[1].getAttributeTag());
	
	EXPECT_EQ(molecules1[2].getAttributeTag(),molecules2[2].getAttributeTag());
	EXPECT_EQ(molecules1[2].getAttributeTag(),molecules3[2].getAttributeTag());
	
	EXPECT_EQ(molecules1[3].getAttributeTag(),molecules2[3].getAttributeTag());
	EXPECT_EQ(molecules1[3].getAttributeTag(),molecules3[3].getAttributeTag());
	
	EXPECT_EQ(molecules1[4].getAttributeTag(),molecules2[4].getAttributeTag());
	EXPECT_EQ(molecules1[4].getAttributeTag(),molecules3[4].getAttributeTag());
	
}

