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
 * @brief Tests for the classes UpdaterReadBfmFile and ReadFullBFMFile
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureJumps.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class ReadBfmFileTest: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureMoleculesIO,FeatureJumps) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> MyIngredients;

  MyIngredients ingredients;

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

TEST_F(ReadBfmFileTest, UpdaterReadBfmFile)
{
  //add a UpdaterReadBfmFile updater to a taskmanager
  TaskManager taskmanager;

  //Constructor of UpdaterReadBfmFile prepares the read-in of HEADER and 1.st frame
  taskmanager.addUpdater(new UpdaterReadBfmFile<MyIngredients>("tests/readbfmfile.test",ingredients,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE));




  //let the taskmanager initalize -> 1.st Frame
  //read-in of HEADER and 1.st frame
  taskmanager.initialize();


  //check if number of monomers was read correctly on initializing the updater
  EXPECT_EQ(10,ingredients.getMolecules().size());
//
  //check some positions
  VectorInt3 position2; position2.setAllCoordinates(5,5,5);
  VectorInt3 position3; position3.setAllCoordinates(3,5,4);
  VectorInt3 position4; position4.setAllCoordinates(5,7,5);
  EXPECT_EQ(position2,ingredients.getMolecules()[2]);
  EXPECT_EQ(position3,ingredients.getMolecules()[3]);
  EXPECT_EQ(position4,ingredients.getMolecules()[4]);

  //check some bond information:

  //number of bonds
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(0));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(3));
  EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(9));
  EXPECT_EQ(1, ingredients.getMolecules().getNumLinks(6));
//   bond partners
  EXPECT_EQ(0, ingredients.getMolecules().getNeighborIdx(6,0));
  EXPECT_EQ(6, ingredients.getMolecules().getNeighborIdx(0,0));
  EXPECT_EQ(1, ingredients.getMolecules().getNeighborIdx(0,1));

  //let the taskmanager run -> 2.nd Frame
  taskmanager.run();

  //now check if the last mcs was read correctly
  EXPECT_EQ(30, ingredients.getMolecules().getAge());
  position2.setAllCoordinates(5+234,5,5-2*67);
  EXPECT_EQ(position2,ingredients.getMolecules()[2]);
  position3.setAllCoordinates(3+234,5,4-2*67);
  EXPECT_EQ(position3,ingredients.getMolecules()[3]);

  //since this updater is supposed to run through the file step by step,
  //it should have executed more than once
  EXPECT_TRUE(1<taskmanager.getNCircles());

}

//this test is the same as the one above, except it used the ReadFullBFMFile
//updater.
TEST_F(ReadBfmFileTest, ReadFullBFMFile)
{
  //add a UpdaterReadBfmFile updater to a taskmanager
  TaskManager taskmanager;
  //Constructor of UpdaterReadBfmFile prepares the read-in of HEADER and last frame
  taskmanager.addUpdater(new UpdaterReadBfmFile<MyIngredients>("tests/readbfmfile.test",ingredients,UpdaterReadBfmFile<MyIngredients>::READ_LAST_CONFIG_SAVE));

  //Constructor of UpdaterReadBfmFile read-in the HEADER and 1.st frame
  //check if number of monomers was read correctly on initializing the updater

//
  //check some positions
  VectorInt3 position2;
  VectorInt3 position3;
  VectorInt3 position4;


  //let the taskmanager initalize -> last Frame
  //read-in of HEADER and last frame
   taskmanager.initialize();

   //now check if the last mcs was read correctly
   EXPECT_EQ(30, ingredients.getMolecules().getAge());
   position2.setAllCoordinates(5+234,5,5-2*67);
   EXPECT_EQ(position2,ingredients.getMolecules()[2]);
   position3.setAllCoordinates(3+234,5,4-2*67);
   EXPECT_EQ(position3,ingredients.getMolecules()[3]);

     //check some bond information:

      //number of bonds
      EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(0));
      EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(3));
      EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(9));
    //   bond partners
      EXPECT_EQ(0, ingredients.getMolecules().getNeighborIdx(6,0));
      EXPECT_EQ(6, ingredients.getMolecules().getNeighborIdx(0,0));
      EXPECT_EQ(1, ingredients.getMolecules().getNeighborIdx(0,1));



  //let the taskmanager run -> last Frame
  taskmanager.run();

  //now check if the last mcs was read correctly
  EXPECT_EQ(30, ingredients.getMolecules().getAge());
  position2.setAllCoordinates(5+234,5,5-2*67);
  EXPECT_EQ(position2,ingredients.getMolecules()[2]);
  position3.setAllCoordinates(3+234,5,4-2*67);
  EXPECT_EQ(position3,ingredients.getMolecules()[3]);

  //since this updater is supposed to run through the complete file on execution
  //it should have executed exactly once.
  EXPECT_EQ(1,taskmanager.getNCircles());

}

TEST_F(ReadBfmFileTest, getMinMaxAge)
{
    //Import the file using the class FileImport
    FileImport<MyIngredients> file("tests/fileImportTest3.test",ingredients);

    //scan file for !mcs and read-in first frame
    file.initialize();

    EXPECT_EQ(2000000,file.getMaxAge());
    EXPECT_EQ(1000,file.getMinAge());
}
