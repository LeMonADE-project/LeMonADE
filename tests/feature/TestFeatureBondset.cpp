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
 * @brief Tests for the bondset featur
 *
 * @date 30.06.2014
 * */
/*****************************************************************************/

#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>

#include <LeMonADE/utility/FastBondset.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/GenerateMonomerType.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class FeatureBondsetTest: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureBox) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;

  //dummy move class used to check response to unknown move type
  class UnknownMove:public MoveBase
  {
	  public:
		template<class IngredientsType> bool check(IngredientsType& ingredients) const
		{
		return ingredients.checkMove(ingredients,*this);
		}

		template<class IngredientsType> void apply(IngredientsType& ingredients)
		{
		ingredients.applyMove(ingredients,*this);
		}

		template <class IngredientsType> void init(const IngredientsType& ingredients){};
  };

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

/************************************************************************/
//test the getBondset and modifyBondset methods of class FeatureBondset
/************************************************************************/
TEST_F(FeatureBondsetTest, GetterAndSetter){
  FeatureBondset<> feature;

  //add a bond
  feature.modifyBondset().addBond(2,2,1,77);

  //check if getBondset and modifyBondset return the same bondset
  EXPECT_EQ(feature.getBondset().getBondVector(77),feature.modifyBondset().getBondVector(77));
  EXPECT_EQ(feature.getBondset().getBondIdentifier(2,2,1),77);

}

/************************************************************************/
//test the exportRead method of FeatureBondset
//for this purpose a testfile with the name "featureBondsetTest.test"
//has to be located in tests source-directory
/************************************************************************/
TEST_F(FeatureBondsetTest, ExportRead){

  //create a file import class to export the read functionality to
  //exportRead is automatically called on creating an instance of FileImport
  Ing ingredients;

  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.setBoxX(10);
  ingredients.setBoxY(10);
  ingredients.setBoxZ(10);

  FileImport<Ing> file ("tests/featureBondsetTest.test",ingredients);

  //scan file for !mcs and read-in first frame
  file.initialize();

  //test if values are read correctly from file
  //these are the vectors defined in the testfile featureBondsetTest.test

  VectorInt3 bondvector1(2,2,1); int32_t identifier1=77;
  VectorInt3 bondvector2(-2,0,-1); int32_t identifier2=45;
  VectorInt3 bondvector3(3,0,-1); int32_t identifier3=117;

  EXPECT_EQ(ingredients.getBondset().getBondVector(77),bondvector1);
  EXPECT_EQ(ingredients.getBondset().getBondVector(45),bondvector2);
  EXPECT_EQ(ingredients.getBondset().getBondVector(117),bondvector3);

  EXPECT_EQ(ingredients.getBondset().getBondIdentifier(2,2,1),77);
  EXPECT_EQ(ingredients.getBondset().getBondIdentifier(-2,0,-1),45);
  EXPECT_EQ(ingredients.getBondset().getBondIdentifier(3,0,-1),117);

}

/***********************************************************************/
//checkMove for unknown moves should return true
/***********************************************************************/
TEST_F(FeatureBondsetTest,checkUnknownMove)
{
	UnknownMove someMove;
	FeatureBondset<> feature;
	EXPECT_TRUE(someMove.check(feature));
}

/***********************************************************************/
//checkMove for MoveLocalBase
/***********************************************************************/
TEST_F(FeatureBondsetTest,checkLocalScBfmMove)
{
	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(10);
	ingredients.setBoxY(10);
	ingredients.setBoxZ(10);

	ingredients.modifyMolecules().resize(2);
	ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
	ingredients.modifyMolecules()[1].setAllCoordinates(2,0,0);
	ingredients.modifyMolecules().connect(0,1);

	ingredients.modifyBondset().addBond(2,0,0,17);
	ingredients.modifyBondset().addBond(-2,0,0,20);

	ingredients.synchronize(ingredients);

	MoveLocalSc move;
	move.init(ingredients);

	//as of now there should be no allowed bonds->any move should be rejected
	//check this for 10 arbitrarily drawn moves.
	for(int i=0;i<10;++i)
	{
		move.init(ingredients);
		EXPECT_FALSE(move.check(ingredients));
	}

	//now disconnect the two monomers and all moves should be accepted
	ingredients.modifyMolecules().disconnect(0,1);
	for(int i=0;i<10;++i)
	{
		move.init(ingredients);
		EXPECT_TRUE(move.check(ingredients));
	}

	//remove the list of bond-vectors
	ingredients.modifyBondset().clear();

	//now introduce bfm classic bondset and check some specific moves
	ingredients.modifyBondset().addBFMclassicBondset();
	ingredients.synchronize(ingredients);

	ingredients.modifyMolecules().connect(0,1);


	for(int i=0;i<10;i++)
	{

		//move of monomer 0 in positive in x-direction must be forbidden
		while( move.getDir().getX()!=1 || move.getIndex()!=0)
		{
			move.init(ingredients);

		}
 		EXPECT_FALSE(move.check(ingredients));

		//move of monomer 1 in negative in x-direction must be forbidden
		while( move.getDir().getX()!=-1 || move.getIndex()!=1)
		{
			move.init(ingredients);

		}
 		EXPECT_FALSE(move.check(ingredients));

		//any move in y and z direction is allowed
		while(move.getDir().getY()==0)
		{
			move.init(ingredients);
		}
		EXPECT_TRUE(move.check(ingredients));

		//any move in y and z direction is allowed
		while(move.getDir().getZ()==0)
		{
			move.init(ingredients);
		}
		EXPECT_TRUE(move.check(ingredients));
	}
}

/************************************************************************/
//test the class ReadBondset
/************************************************************************/
TEST_F(FeatureBondsetTest, ReadBondset){


  //define some input streams for the ReadBondset to read from
  istringstream inputStandard;
  istringstream inputSpaces;
  istringstream inputWrongFormat;
  istringstream inputMissingComponent;
  istringstream inputMissingIdentifier;

  //newline at beginning necessary, because the command structure is
  // !set_of_bondvectors
  // x y z:id
  //and the parser leaves the get pointer behind the initial keyword
  inputStandard.str("\n2 2 1:77");
  inputSpaces.str("  \n  2    1   2  :    78");
  inputWrongFormat.str("\n1 2 2-79");
  inputMissingIdentifier.str("\n-2 2 1 80");
  inputMissingComponent.str("\n2 1:81");

  //create a reference bondvector
  VectorInt3 bondvector1(2,2,1);
  int32_t identifier1=77;

  //create and test the command object
  FeatureBondset<> bondset;
  ReadBondset<FeatureBondset<> > command(bondset);

  //check reaction to correct standard input
  command.setInputStream(&inputStandard);
  command.execute();
  EXPECT_EQ(bondset.getBondset().getBondVector(77),bondvector1);
  EXPECT_EQ(bondset.getBondset().getBondIdentifier(2,2,1),77);

  //check reaction to input containing extra whitespaces
  command.setInputStream(&inputSpaces);
  command.execute();
  VectorInt3 bondvector2(2,1,2);
  int32_t identifier2=78;
  EXPECT_EQ(bondset.getBondset().getBondVector(78),bondvector2);
  EXPECT_EQ(bondset.getBondset().getBondIdentifier(2,1,2),78);

  //check i reaction to input with wrong format of different types
  command.setInputStream(&inputWrongFormat);
  EXPECT_THROW(command.execute(),std::runtime_error);

  command.setInputStream(&inputMissingComponent);
  EXPECT_THROW(command.execute(),std::runtime_error);

  command.setInputStream(&inputMissingIdentifier);
  EXPECT_THROW(command.execute(),std::runtime_error);

}

/************************************************************************/
//test the class ReadBondset
/************************************************************************/
TEST_F(FeatureBondsetTest, WriteBondsetTest){

  //the test sets up a bondset, writes it to a file,
  //reads the file back in to a second bondset, and compares the two sets.
  Ing ingredients;
  Ing bondsetCheck;

  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.setBoxX(10);
  ingredients.setBoxY(10);
  ingredients.setBoxZ(10);

  bondsetCheck.setPeriodicX(true);
  bondsetCheck.setPeriodicY(true);
  bondsetCheck.setPeriodicZ(true);
  bondsetCheck.setBoxX(10);
  bondsetCheck.setBoxY(10);
  bondsetCheck.setBoxZ(10);
  FileImport<Ing> file ("tests/featureBondsetTest.test",ingredients);

  //scan file for !mcs and read-in first frame
  file.initialize();

  //add some bonds
  ingredients.modifyBondset().addBond(1,1,1,101);
  ingredients.modifyBondset().addBond(2,2,2,102);
  ingredients.modifyBondset().addBond(3,3,3,103);

  //create output file
  WriteBondset<Ing> output(ingredients);
  ofstream outfile;
  outfile.open("tmpBondSet.bfm",ios_base::app|ios_base::binary);
  output.writeStream(outfile);
  outfile.close();

  //read file back in
  FileImport<Ing> infile("tmpBondSet.bfm",bondsetCheck);

  //scan file for !mcs and read-in first frame
  infile.initialize();


  //compare the two sets
  VectorInt3 bond1(1,1,1);
  int32_t identifier1=101;
  VectorInt3 bond2(2,2,2);
  int32_t identifier2=102;
  VectorInt3 bond3(3,3,3);
  int32_t identifier3=103;

  EXPECT_EQ(bond1,(bondsetCheck.getBondset().getBondVector(identifier1)));
  //EXPECT_EQ(101,(bondsetCheck.getBondset().begin()->GetW()));

  EXPECT_EQ(bond2,((bondsetCheck.getBondset().getBondVector(identifier2))));
 // EXPECT_EQ(102,(++(bondsetCheck.getBondset().begin()))->GetW());

  EXPECT_EQ(bond3,((bondsetCheck.getBondset().getBondVector(identifier3))));
  //EXPECT_EQ(103,(--(bondsetCheck.getBondset().end()))->GetW());

  //remove the temporary file
  EXPECT_EQ(0,remove("tmpBondSet.bfm"));
}


