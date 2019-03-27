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
 * @brief Tests for the classes AnalyzerWriteBfmFile
 *
 * @author Martin
 * @date 18.06.2014
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/utility/Vector3D.h>


using namespace std;
// class with body that is needed for every test, only used by TEST_F()
/*****************************************************************************/
/**
 * @class WriteBfmFileTest
 * @brief prepare system and check functions for AnalyzerWriteBfmFile Analyzer Test
 * @details when using TEST_F(WriteBfmFileTest, -testname-), an instance of
 * Ingredients (ingredients) with features Box, Bondset and Attributes is prepared.
 * To check constructor, GoogleTest needs void functions, that are also provided
 * in this class \n
 * for large/less output just comment/uncomment the depending parts of the test classes
 * @todo: check FeatureJumps ()?, FeatureSolvent and !add_bonds command
 * */
/*****************************************************************************/
class WriteBfmFileTest: public ::testing::Test{
protected:
  //void functions to check constructor(not possible to make explicit contructur assertion
  void checkflags(string _filename){
    AnalyzerWriteBfmFile<MyIngredients>BfmWriter(_filename, ingredients,4);
    BfmWriter.initialize();
  }
  void checkexistingfile(string _filename){
    AnalyzerWriteBfmFile<MyIngredients>BfmWriter(_filename, ingredients, AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
    BfmWriter.initialize();
  }
  void checkfileopener(string _filename){
    AnalyzerWriteBfmFile<MyIngredients>BfmWriter(_filename, ingredients, AnalyzerWriteBfmFile<MyIngredients>::APPEND);
    BfmWriter.initialize();
  }

  //define system
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients < Config> MyIngredients;
  MyIngredients ingredients;

  /* suppress cout output for better readability -->un/comment here:*/
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
/**
 * @fn TEST_F(WriteBfmFileTest, CheckConstructor)
 * @brief Test function for constructor of AnalyzerWriteBfmFile analyzer.
 * */
/*****************************************************************************/
TEST_F(WriteBfmFileTest, CheckConstructor)
{
  /*
  //create a system using FeatureBondset and the Molecules template
  typedef LOKI_TYPELIST_3(FeatureBox,FeatureBondset,FeatureAttributes) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients < Config> MyIngredients;
  MyIngredients ingredients;
  */
  ingredients.modifyMolecules().addMonomer(42,13,40);
  //use FeatureBox
  ingredients.setBoxX(128);
  ingredients.setBoxY(128);
  ingredients.setBoxZ(128);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(false);

  ingredients.synchronize(ingredients);

  //produce some output
  string filename("tests/existing_file.test");
  AnalyzerWriteBfmFile<MyIngredients> BfmWriter(filename, ingredients);
  BfmWriter.initialize();
  BfmWriter.execute();
  BfmWriter.cleanup();

  //check if Flags are read correctly (eigther APPEND=0 or NEWFILE=1)
  EXPECT_THROW(checkflags(filename), std::runtime_error);

  // check if writer finds existing files correctly
  EXPECT_THROW(checkexistingfile(filename), std::runtime_error);

  // check fileopener
  EXPECT_NO_THROW(checkfileopener(filename));

  remove(filename.c_str());
}

/*****************************************************************************/
/**
 * @fn TEST_F(WriteBfmFileTest, CheckOutput)
 * @brief Test function for correct output of AnalyzerWriteBfmFile analyzer.
 * */
/*****************************************************************************/
TEST_F(WriteBfmFileTest, CheckOutput)
{
  //use !mcs
  ingredients.modifyMolecules().addMonomer(42,13,40);
  ingredients.modifyMolecules().addMonomer(41,15,40);
  ingredients.modifyMolecules().addMonomer(40,17,39);
  ingredients.modifyMolecules().addMonomer(45,13,41);
  //use !bondset
  ingredients.modifyBondset().addBond(-1,2,0,38);
  ingredients.modifyBondset().addBond(-1,2,-1,67);
  ingredients.modifyBondset().addBond(3,0,1,114);

  ingredients.modifyBondset().addBond(1,-2,0,41);
  ingredients.modifyBondset().addBond(1,-2,1,52);
  ingredients.modifyBondset().addBond(-3,0,-1,123);
  //use !bonds
  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().connect(0,3);
  //use !mcs=age
  ingredients.modifyMolecules().setAge(5854);
  //use !box_x _y _z
  ingredients.setBoxX(128);
  ingredients.setBoxY(256);
  ingredients.setBoxZ(64);
  //use !periodic_x _y _z
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(false);
  //use attributes
  ingredients.modifyMolecules()[0].setAttributeTag(11);
  ingredients.modifyMolecules()[1].setAttributeTag(11);
  ingredients.modifyMolecules()[2].setAttributeTag(11);
  ingredients.modifyMolecules()[3].setAttributeTag(22);

  ingredients.synchronize(ingredients);

  //prepare output file
  string filename("tests/writebfmfile.test");

  /* using explicit definition and execution of writer and reader.
   * taskmanager not possible: wrong order of writer and reader  */
  AnalyzerWriteBfmFile<MyIngredients> BfmWriter(filename, ingredients, AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
  BfmWriter.initialize();
  BfmWriter.execute();
  BfmWriter.cleanup();

  //check file by read it in
  MyIngredients iningredients;
  UpdaterReadBfmFile<MyIngredients> BfmReader(filename, iningredients,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE);
  BfmReader.initialize();
  BfmReader.execute();
  BfmReader.cleanup();
  //check !number_of_monomers
  EXPECT_EQ(4,iningredients.getMolecules().size());
  //check age
  EXPECT_EQ(5854,iningredients.getMolecules().getAge());
  //check box
  EXPECT_EQ(128,iningredients.getBoxX());
  EXPECT_EQ(256,iningredients.getBoxY());
  EXPECT_EQ(64,iningredients.getBoxZ());
  //check periodocity
  EXPECT_TRUE(iningredients.isPeriodicX());
  EXPECT_TRUE(iningredients.isPeriodicY());
  EXPECT_FALSE(iningredients.isPeriodicZ());
  //check connections
  EXPECT_TRUE(iningredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(iningredients.getMolecules().areConnected(1,2));
  EXPECT_TRUE(iningredients.getMolecules().areConnected(0,3));
  EXPECT_FALSE(iningredients.getMolecules().areConnected(1,3));
  //check some monomer coordinates
  EXPECT_EQ(42,iningredients.getMolecules()[0].getX());
  EXPECT_EQ(15,iningredients.getMolecules()[1].getY());
  EXPECT_EQ(39,iningredients.getMolecules()[2].getZ());
  EXPECT_EQ(45,iningredients.getMolecules()[3].getX());
  //check some bondset coordinates
  EXPECT_EQ(-1,iningredients.getBondset().getBondVector(38).getX());
  EXPECT_EQ(-1,iningredients.getBondset().getBondVector(67).getZ());
  EXPECT_EQ(3,iningredients.getBondset().getBondVector(114).getX());
  //check attributes
  EXPECT_EQ(11,iningredients.getMolecules()[0].getAttributeTag());
  EXPECT_EQ(11,iningredients.getMolecules()[1].getAttributeTag());
  EXPECT_EQ(11,iningredients.getMolecules()[2].getAttributeTag());
  EXPECT_EQ(22,iningredients.getMolecules()[3].getAttributeTag());

  remove(filename.c_str());
}
