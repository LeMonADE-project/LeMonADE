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

/**
 * @file
 * @brief Tests for the feature FeatureMoleculesIO
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <stdexcept>
#include <cstdio>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>

using namespace std;

/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class FeatureMoleculesIOTest: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_1(FeatureMoleculesIO) Features;
  typedef ConfigureSystem<VectorInt3,Features,4> Config;
  typedef Ingredients < Config> MyIngredients;

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


//test the if the read commands !number_of_monomers, !bonds and !mcs
//are exported correctly and work.
TEST_F(FeatureMoleculesIOTest, ExportRead){
  MyIngredients ingredients;

  //read a testfile
  FileImport<MyIngredients> file("tests/molecules.test",ingredients);

  //check if number of monomers was read correctly
  file.read();
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
//   bond partners
  EXPECT_EQ(0, ingredients.getMolecules().getNeighborIdx(5,0));
  EXPECT_EQ(5, ingredients.getMolecules().getNeighborIdx(0,0));
  EXPECT_EQ(1, ingredients.getMolecules().getNeighborIdx(0,1));

  //check age (first mcs)
  EXPECT_EQ(10, ingredients.getMolecules().getAge());

  //check age (next mcs)
  file.read();
  EXPECT_EQ(20, ingredients.getMolecules().getAge());

  //in the next mcs there are jumps, which is not part of the standard.
  //an exception should be thrown
  EXPECT_THROW(file.read(),std::runtime_error);
  /* **** this part of the test was before FeatureJumps was introduced, i.e. when
   * jumps were still part of the standart !mcs command. it is left here commented
   * out because it may still be used in a test for FeatureJump
  //check age (next mcs)
  file.read();
  EXPECT_EQ(30, ingredients.getMolecules().getAge());
  position2.setAllCoordinates(5+234,5,5-2*67);
  EXPECT_EQ(position2,ingredients.getMolecules()[2]);
  position3.setAllCoordinates(3+234,5,4-2*67);
  EXPECT_EQ(position3,ingredients.getMolecules()[3]);

  //check age (next mcs)
  EXPECT_THROW(file.read(),std::runtime_error);
  EXPECT_EQ(40, ingredients.getMolecules().getAge());
*/

}

//test if the read commands !number_of_monomers, !bondmovecounters and !mcs
//are exported correctly and work
TEST_F(FeatureMoleculesIOTest, ExportWrite){
  MyIngredients ingredients;

  //set some values on ingredients
  //ingredients.modifyMolecules().resize(9);
  ingredients.modifyMolecules().setAge(10000);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.setBoxX(128);
  ingredients.setBoxY(128);
  ingredients.setBoxZ(128);

  VectorInt3 monomer1; monomer1.setAllCoordinates(1,1,1);
  VectorInt3 monomer2; monomer2.setAllCoordinates(2,2,2);
  VectorInt3 monomer3; monomer3.setAllCoordinates(3,3,3);
  VectorInt3 monomer4; monomer4.setAllCoordinates(1,4,1);
  VectorInt3 monomer5; monomer5.setAllCoordinates(2,5,2);
  VectorInt3 monomer6; monomer6.setAllCoordinates(3,6,3);
  VectorInt3 monomer7; monomer7.setAllCoordinates(7,7,7);
  VectorInt3 monomer8; monomer8.setAllCoordinates(6,9,8);
  VectorInt3 monomer9; monomer9.setAllCoordinates(6,12,8);
  VectorInt3 monomer10; monomer10.setAllCoordinates(7,13,9);

  ingredients.modifyMolecules().addMonomer(monomer1);
  ingredients.modifyMolecules().addMonomer(monomer2);
  ingredients.modifyMolecules().addMonomer(monomer3);
  ingredients.modifyMolecules().addMonomer(monomer4);
  ingredients.modifyMolecules().addMonomer(monomer5);
  ingredients.modifyMolecules().addMonomer(monomer6);
  ingredients.modifyMolecules().addMonomer(monomer7);
  ingredients.modifyMolecules().addMonomer(monomer8);
  ingredients.modifyMolecules().addMonomer(monomer9);
  ingredients.modifyMolecules().addMonomer(monomer10);

  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().connect(3,4);
  ingredients.modifyMolecules().connect(4,5);
  ingredients.modifyMolecules().connect(6,7);
  ingredients.modifyMolecules().connect(7,8);
  ingredients.modifyMolecules().connect(8,9);
  ingredients.modifyMolecules().connect(1,4);


  ingredients.modifyBondset().addBond(1,1,1,')');
  ingredients.modifyBondset().addBond(0,3,0,'o');
  ingredients.modifyBondset().addBond(-1,2,1,':');

  //now create file with the information entered above
  string filename("tmpMoleculesRW.bfm");
  AnalyzerWriteBfmFile<MyIngredients> outputFile(filename,ingredients,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
  outputFile.initialize();
  outputFile.execute();
  outputFile.closeFile();

  cout<<"start reading"<<endl;
  //now read the file back in and compare
  MyIngredients checkIngredients;
  FileImport<MyIngredients> inputFile(filename,checkIngredients);
  inputFile.initialize();

  //check age
  EXPECT_EQ(checkIngredients.getMolecules().getAge(),10000);
  //check number of monomers
  EXPECT_EQ(checkIngredients.getMolecules().size(),10);
  //check positions of monomers
  EXPECT_EQ(monomer1,checkIngredients.getMolecules()[0]);
  EXPECT_EQ(monomer9,checkIngredients.getMolecules()[8]);
  //check connections
  EXPECT_TRUE(checkIngredients.getMolecules().areConnected(0,1));
  EXPECT_TRUE(checkIngredients.getMolecules().areConnected(4,5));
  EXPECT_TRUE(checkIngredients.getMolecules().areConnected(4,1));
  EXPECT_FALSE(checkIngredients.getMolecules().areConnected(2,3));
  //remove temporary file tmp.bfm (and let test fail if not removed)
  EXPECT_EQ(0,remove(filename.c_str()));

}



TEST_F(FeatureMoleculesIOTest,CompressedSolventLowDensity)
{
	MyIngredients ingredients;

	//set some values on ingredients
	ingredients.modifyMolecules().setAge(10000);
	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(32);
	ingredients.setBoxY(32);
	ingredients.setBoxZ(64);

	//set 500 solvents in the box
	int offset=32;
	std::map<int32_t,int32_t> coordinates;
	for(int n=0;n<500;n++)
	{
		int x=offset+n;
		int y=offset+2*n;
		int z=offset+3*n;
		VectorInt3 monomer1; monomer1.setAllCoordinates(x,y,z);
		ingredients.modifyMolecules().addMonomer(monomer1);
		//save folded positions in a map
		//transform the position of the solvent to a linear index
		int32_t foldedX=x%32;
		int32_t foldedY=y%32;
		int32_t foldedZ=z%64;
		int32_t linIndex=foldedZ*32*32+foldedY*32+foldedX;
		coordinates[linIndex]+=1;

	}

	ingredients.setCompressedOutputIndices(0,499);
	ingredients.synchronize(ingredients);

	//now create file with the information entered above
	string filename("tmpMoleculesSoventRW.bfm");
	AnalyzerWriteBfmFile<MyIngredients> outputFile(filename,ingredients,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
	outputFile.initialize();
	outputFile.execute();
	outputFile.closeFile();


	MyIngredients testIngredients;
	UpdaterReadBfmFile<MyIngredients> inputFile(filename,testIngredients,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE);
	inputFile.initialize();
	inputFile.execute();
	inputFile.closeFile();

	//construct a map with folded coordinates again
	std::map<int32_t,int32_t> testCoordinates;
	for(size_t i=0;i<testIngredients.getMolecules().size();i++)
	{
		//transform the position of the solvent to a linear index
		int32_t foldedX=testIngredients.getMolecules()[i].getX()%32;
		int32_t foldedY=testIngredients.getMolecules()[i].getY()%32;
		int32_t foldedZ=testIngredients.getMolecules()[i].getZ()%64;
		int32_t linIndex=foldedZ*32*32+foldedY*32+foldedX;
		testCoordinates[linIndex]+=1;
	}
	//sovents should be same range in both
	EXPECT_EQ(ingredients.getCompressedOutputIndices(),testIngredients.getCompressedOutputIndices());
	//maps with folded indices should also be equal
	EXPECT_EQ(coordinates,testCoordinates);

	//remove temporary file  (and let test fail if not removed)
	EXPECT_EQ(0,remove(filename.c_str()));

}

TEST_F(FeatureMoleculesIOTest,CompressedSolventLowDensitySplit)
{
	MyIngredients ingredients;

	//set some values on ingredients
	ingredients.modifyMolecules().setAge(10000);
	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(32);
	ingredients.setBoxY(64);
	ingredients.setBoxZ(64);

	//set 500 solvents in the box
	int offset=32;
	std::map<int32_t,int32_t> coordinates;
	for(int n=0;n<500;n++)
	{
		int x=offset+n;
		int y=offset+2*n;
		int z=offset+3*n;
		VectorInt3 monomer1; monomer1.setAllCoordinates(x,y,z);
		ingredients.modifyMolecules().addMonomer(monomer1);
		//save folded positions in a map
		//transform the position of the solvent to a linear index
		int32_t foldedX=x%32;
		int32_t foldedY=y%64;
		int32_t foldedZ=z%64;
		int32_t linIndex=foldedZ*32*64+foldedY*32+foldedX;
		coordinates[linIndex]+=1;

	}

	ingredients.setCompressedOutputIndices(10,20);
	ingredients.setCompressedOutputIndices(30,480);
	ingredients.synchronize(ingredients);

	//now create file with the information entered above
	string filename("tmpMoleculesSoventRW.bfm");
	AnalyzerWriteBfmFile<MyIngredients> outputFile(filename,ingredients,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
	outputFile.initialize();
	outputFile.execute();
	outputFile.closeFile();


	MyIngredients testIngredients;
	UpdaterReadBfmFile<MyIngredients> inputFile(filename,testIngredients,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE);
	inputFile.initialize();
	inputFile.execute();
	inputFile.closeFile();

	//construct a map with folded coordinates again
	std::map<int32_t,int32_t> testCoordinates;
	for(size_t i=0;i<testIngredients.getMolecules().size();i++)
	{
		//transform the position of the solvent to a linear index
		int32_t foldedX=testIngredients.getMolecules()[i].getX()%32;
		int32_t foldedY=testIngredients.getMolecules()[i].getY()%64;
		int32_t foldedZ=testIngredients.getMolecules()[i].getZ()%64;
		int32_t linIndex=foldedZ*32*64+foldedY*32+foldedX;
		testCoordinates[linIndex]+=1;
	}
	//sovents should be same range in both
	EXPECT_EQ(ingredients.getCompressedOutputIndices(),testIngredients.getCompressedOutputIndices());
	//maps with folded indices should also be equal
	EXPECT_EQ(coordinates,testCoordinates);

	//remove temporary file  (and let test fail if not removed)
	EXPECT_EQ(0,remove(filename.c_str()));

}



TEST_F(FeatureMoleculesIOTest,CompressedSolventHighDensity)
{
	MyIngredients ingredients;

	//set some values on ingredients
	ingredients.modifyMolecules().setAge(10000);
	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(32);
	ingredients.setBoxY(32);
	ingredients.setBoxZ(64);

	//set 500 solvents in the box
	int offset=32;
	std::map<int32_t,int32_t> coordinates;
	for(int n=0;n<2000;n++)
	{
		int x=offset+n;
		int y=offset+2*n;
		int z=offset+3*n;
		VectorInt3 monomer1; monomer1.setAllCoordinates(x,y,z);
		ingredients.modifyMolecules().addMonomer(monomer1);
		//save folded positions in a map
		//transform the position of the solvent to a linear index
		int32_t foldedX=x%32;
		int32_t foldedY=y%32;
		int32_t foldedZ=z%64;
		int32_t linIndex=foldedZ*32*32+foldedY*32+foldedX;
		coordinates[linIndex]+=1;

	}

	ingredients.setCompressedOutputIndices(0,1999);
	ingredients.synchronize(ingredients);

	//now create file with the information entered above
	string filename("tmpMoleculesSoventRW.bfm");
	AnalyzerWriteBfmFile<MyIngredients> outputFile(filename,ingredients,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
	outputFile.initialize();
	outputFile.execute();
	outputFile.closeFile();


	MyIngredients testIngredients;
	UpdaterReadBfmFile<MyIngredients> inputFile(filename,testIngredients,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE);
	inputFile.initialize();
	inputFile.execute();
	inputFile.closeFile();

	//construct a map with folded coordinates again
	std::map<int32_t,int32_t> testCoordinates;
	for(size_t i=0;i<testIngredients.getMolecules().size();i++)
	{
		//transform the position of the solvent to a linear index
		int32_t foldedX=testIngredients.getMolecules()[i].getX()%32;
		int32_t foldedY=testIngredients.getMolecules()[i].getY()%32;
		int32_t foldedZ=testIngredients.getMolecules()[i].getZ()%64;
		int32_t linIndex=foldedZ*32*32+foldedY*32+foldedX;
		testCoordinates[linIndex]+=1;
	}
	//sovents should be same range in both
	EXPECT_EQ(ingredients.getCompressedOutputIndices(),testIngredients.getCompressedOutputIndices());
	//maps with folded indices should also be equal
	EXPECT_EQ(coordinates,testCoordinates);

	//remove temporary file  (and let test fail if not removed)
	EXPECT_EQ(0,remove(filename.c_str()));

}

TEST_F(FeatureMoleculesIOTest,CompressedSolventHighDensitySplit)
{
	MyIngredients ingredients;

	//set some values on ingredients
	ingredients.modifyMolecules().setAge(10000);
	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(32);
	ingredients.setBoxY(64);
	ingredients.setBoxZ(64);

	//set 500 solvents in the box
	int offset=32;
	std::map<int32_t,int32_t> coordinates;
	for(int n=0;n<2000;n++)
	{
		int x=offset+n;
		int y=offset+2*n;
		int z=offset+3*n;
		VectorInt3 monomer1; monomer1.setAllCoordinates(x,y,z);
		ingredients.modifyMolecules().addMonomer(monomer1);
		//save folded positions in a map
		//transform the position of the solvent to a linear index
		int32_t foldedX=x%32;
		int32_t foldedY=y%64;
		int32_t foldedZ=z%64;
		int32_t linIndex=foldedZ*32*64+foldedY*32+foldedX;
		coordinates[linIndex]+=1;

	}

	ingredients.setCompressedOutputIndices(0,20);
	ingredients.setCompressedOutputIndices(30,280);
	ingredients.setCompressedOutputIndices(300,1999);
	ingredients.synchronize(ingredients);

	//now create file with the information entered above
	string filename("tmpMoleculesSoventRW.bfm");
	AnalyzerWriteBfmFile<MyIngredients> outputFile(filename,ingredients,AnalyzerWriteBfmFile<MyIngredients>::NEWFILE);
	outputFile.initialize();
	outputFile.execute();
	outputFile.closeFile();


	MyIngredients testIngredients;
	UpdaterReadBfmFile<MyIngredients> inputFile(filename,testIngredients,UpdaterReadBfmFile<MyIngredients>::READ_STEPWISE);
	inputFile.initialize();
	inputFile.execute();
	inputFile.closeFile();

	//construct a map with folded coordinates again
	std::map<int32_t,int32_t> testCoordinates;
	for(size_t i=0;i<testIngredients.getMolecules().size();i++)
	{
		//transform the position of the solvent to a linear index
		int32_t foldedX=testIngredients.getMolecules()[i].getX()%32;
		int32_t foldedY=testIngredients.getMolecules()[i].getY()%64;
		int32_t foldedZ=testIngredients.getMolecules()[i].getZ()%64;
		int32_t linIndex=foldedZ*32*64+foldedY*32+foldedX;
		testCoordinates[linIndex]+=1;
	}
	//sovents should be same range in both
	EXPECT_EQ(ingredients.getCompressedOutputIndices(),testIngredients.getCompressedOutputIndices());
	//maps with folded indices should also be equal
	EXPECT_EQ(coordinates,testCoordinates);

	//remove temporary file  (and let test fail if not removed)
	EXPECT_EQ(0,remove(filename.c_str()));

}


TEST_F(FeatureMoleculesIOTest,getSetCompressedOutputIndices)
{
	MyIngredients ingredients;

	EXPECT_EQ(ingredients.getCompressedOutputIndices().size(),0);

	ingredients.setCompressedOutputIndices(10,112);
	ingredients.setCompressedOutputIndices(200,400);
	EXPECT_EQ(ingredients.getCompressedOutputIndices().size(),2);

	ingredients.setCompressedOutputIndices(10,114);
	EXPECT_EQ(ingredients.getCompressedOutputIndices().at(10),114);
}
