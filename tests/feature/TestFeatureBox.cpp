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

#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveAddMonomerSc.h>

using namespace std;

/*****************************************************************************/
/**
 * @file
 * @brief all tests relevant for FeatureBox
 * @details run this tests only: ./LemonadeTests --gtest_filter=*FeatureBox* \n
 * Includes tests for the class FeatureBox, as well as tests for the
 * classes ReadBox_ and ReadPeriodic_
 * */
/*************************************************************************** */

/*****************************************************************************/
/**
 * @class FeatureBoxTest
 * @brief prepare ingredients system
 * @details when using TEST_F(FeatureBoxTest, -testname-), an instance of
 * Ingredients (ingredients) with features Box and Bondset is prepared. \n
 * for large/less output just comment/uncomment the depending parts of the test classes
 * */
/*************************************************************************** */
class FeatureBoxTest: public ::testing::Test{
protected:
  typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureBox) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;
  Ing ingredients;

public:

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

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, GetterDeath)
 * @brief test for exit of program if box sizes are not set
/*************************************************************************** */
TEST_F(FeatureBoxTest, GetterDeath){
  EXPECT_THROW(ingredients.getBoxX(),std::runtime_error );
  EXPECT_THROW(ingredients.getBoxY(),std::runtime_error);
  EXPECT_THROW(ingredients.getBoxZ(),std::runtime_error);

  /* why no exceptions thrown by assertPeriodicitySet???*/
  EXPECT_THROW(ingredients.isPeriodicX(),std::runtime_error );
  EXPECT_THROW(ingredients.isPeriodicY(),std::runtime_error);
  EXPECT_THROW(ingredients.isPeriodicZ(),std::runtime_error);
  /* */
}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, GetterAndSetter)
 * @brief test getter and setter functions
/*************************************************************************** */
TEST_F(FeatureBoxTest, GetterAndSetter) {
  FeatureBox box;

  //change the streambuffer of cin to mimic user input
  //this is needed for the test because the function
  //error::handleWarning expects user input via cin
//   std::streambuf* original=std::cin.rdbuf();
//   std::istringstream mimicInput("123");
//   std::cin.rdbuf(mimicInput.rdbuf());

  //now do the tests requiring user input
  //test default values of periodicity
  /*
  EXPECT_TRUE(box.isPeriodicX());
  EXPECT_TRUE(box.isPeriodicY());
  EXPECT_TRUE(box.isPeriodicZ());/* */

  //reset the streambuffer of cin to original value
//   std::cin.rdbuf(original);

  //now test periodicity getter and setter
  box.setPeriodicX(false);
  box.setPeriodicY(false);
  box.setPeriodicZ(false);

  EXPECT_FALSE(box.isPeriodicX());
  EXPECT_FALSE(box.isPeriodicY());
  EXPECT_FALSE(box.isPeriodicZ());

  box.setPeriodicX(true);
  box.setPeriodicY(true);
  box.setPeriodicZ(true);

  EXPECT_TRUE(box.isPeriodicX());
  EXPECT_TRUE(box.isPeriodicY());
  EXPECT_TRUE(box.isPeriodicZ());

  //test boxsize getter and setter
  box.setBoxX(100);
  box.setBoxY(200);
  box.setBoxZ(300);

  EXPECT_EQ(100,box.getBoxX());
  EXPECT_EQ(200,box.getBoxY());
  EXPECT_EQ(300,box.getBoxZ());

  //test number of lattice sites
  EXPECT_EQ(100*200*300,box.getNumberOfLatticeSites());

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, ExportRead)
 * @brief test if the export of the read-from-file funcionality is working. \n
 * For this purpose a testfile with the name \b "tests/featureBoxTest.test"
 * has to be located in tests source-directory
/*************************************************************************** */
TEST_F(FeatureBoxTest, ExportRead){

  //initialize a file import with the testfile and feature
  //exportRead is automatically called on creating an instance of FileImport
  //  FeatureBox feature;
  FileImport<Ing> file("tests/featureBoxTest.test",ingredients);

  //scan file for !mcs and read-in first frame
  file.initialize();

  //test if values are read correctly from file
   EXPECT_EQ(ingredients.getBoxX(),100);
   EXPECT_EQ(ingredients.getBoxY(),200);
   EXPECT_EQ(ingredients.getBoxZ(),300);
   EXPECT_TRUE(ingredients.isPeriodicX());
   EXPECT_TRUE(ingredients.isPeriodicY());
   EXPECT_FALSE(ingredients.isPeriodicZ());

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, Synchronize)
 * @brief test if the synchronize function is working correctly.
/*************************************************************************** */
TEST_F(FeatureBoxTest, Synchronize)
{
	ingredients.modifyMolecules().resize(1);
	ingredients.modifyMolecules()[0].setAllCoordinates(20,20,20);

	//trying to synchronize without box size being set in all directions
	//should cause an exception.
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);
	ingredients.setBoxX(10);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);
	ingredients.setBoxY(10);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);
	ingredients.setBoxZ(10);
	//now all box dimenstions are set, but if we make the box non
	//periodic, there should be an exception thrown because the
	//particles are not inside the box

	ingredients.setPeriodicX(false);
	ingredients.setPeriodicY(false);
	ingredients.setPeriodicZ(false);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);
	ingredients.setPeriodicX(true);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);
	ingredients.setPeriodicY(true);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);
	ingredients.setPeriodicZ(true);
	//now the box is periodic and the system should synchronize
	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	//now, if we set the periodicity back, it should throw an exception again
	ingredients.setPeriodicX(false);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(false);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(false);
	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, checkUnknownMove)
 * @brief test the UnknownMove check
/*************************************************************************** */
//unknown moves that derive from the class Move should be accepted
TEST_F(FeatureBoxTest,checkUnknownMove)
{
	UnknownMove someMove;
	FeatureBox feature;
	//check pure feature
	EXPECT_TRUE(someMove.check(feature));
	//check feature in ingredients
	EXPECT_TRUE(someMove.check(ingredients));
}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, checkLocalBfmMoveMove)
 * @brief tests the MoveLocalBase check
/*************************************************************************** */
TEST_F(FeatureBoxTest, checkLocalBfmMove)
{
	ingredients.setBoxX(10);
	ingredients.setBoxY(10);
	ingredients.setBoxZ(10);

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);

	ingredients.modifyMolecules().resize(1);
	ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);


	// coordinates for non-periodicity 0 <= x,y,z <= BOX-2 (Wall at BOX)
	// coordinates for     periodicity 0 <= x,y,z <= BOX-1 (BOX == 0)
	MoveLocalSc move;
	//first,check move to x<0
	move.init(ingredients);
	while(move.getDir().getX()>-1)
	{
		move.init(ingredients);
	}

	EXPECT_TRUE(move.check(ingredients));
	ingredients.setPeriodicX(false);
	EXPECT_FALSE(move.check(ingredients));


	//now check move to x>boxX-1

	ingredients.setPeriodicX(true);
	ingredients.modifyMolecules()[0].setAllCoordinates(ingredients.getBoxX()-2,0,0);
	while(move.getDir().getX()<1)
	{
		move.init(ingredients);
	}
	EXPECT_TRUE(move.check(ingredients));
	ingredients.setPeriodicX(false);
	EXPECT_FALSE(move.check(ingredients));


	//now y direction, first y<0, then y>boxY-1
	while(move.getDir().getY()>-1)
	{
		move.init(ingredients);
	}

	EXPECT_TRUE(move.check(ingredients));
	ingredients.setPeriodicY(false);
	EXPECT_FALSE(move.check(ingredients));

	ingredients.setPeriodicY(true);
	ingredients.modifyMolecules()[0].setAllCoordinates(ingredients.getBoxX()-2,ingredients.getBoxY()-2,0);
	while(move.getDir().getY()<1)
	{
		move.init(ingredients);
	}
	EXPECT_TRUE(move.check(ingredients));
	ingredients.setPeriodicY(false);
	EXPECT_FALSE(move.check(ingredients));

	//now z direction, first z<0, then z>boxZ-1
	while(move.getDir().getZ()>-1)
	{
		move.init(ingredients);
	}

	EXPECT_TRUE(move.check(ingredients));
	ingredients.setPeriodicZ(false);
	EXPECT_FALSE(move.check(ingredients));

	ingredients.setPeriodicZ(true);
	ingredients.modifyMolecules()[0].setAllCoordinates(ingredients.getBoxX()-2,ingredients.getBoxY()-2,ingredients.getBoxZ()-2);
	while(move.getDir().getZ()<1)
	{
		move.init(ingredients);
	}
	EXPECT_TRUE(move.check(ingredients));
	ingredients.setPeriodicZ(false);
	EXPECT_FALSE(move.check(ingredients));
}


/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxTest, checkLocalBfmMoveMove)
 * @brief tests the MoveLocalBase check
/*************************************************************************** */
TEST_F(FeatureBoxTest, checkAddMonomerMove)
{
	ingredients.setBoxX(10);
	ingredients.setBoxY(10);
	ingredients.setBoxZ(10);

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);


	// coordinates for non-periodicity 0 <= x,y,z <= BOX-2 (Wall at BOX)

	// first check with periodic boundary conditions
	MoveAddMonomerSc<> move;

	move.init(ingredients);
	move.setPosition(-1,3,5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,-3,5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,3,-5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(13,-3,5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,30,5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(2,3,50);
	EXPECT_TRUE(move.check(ingredients));

	//now unset periodic boundary conditions in x direction
	ingredients.setPeriodicX(false);

	move.init(ingredients);
	move.setPosition(-1,3,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(0,3,5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(15,3,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(10,3,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(9,3,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(8,3,5);
	EXPECT_TRUE(move.check(ingredients));

	//now unset periodic boundary conditions  in y direction
	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(false);

	move.init(ingredients);
	move.setPosition(1,-3,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(0,0,5);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,13,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,10,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,9,5);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,8,5);
	EXPECT_TRUE(move.check(ingredients));

	//now unset periodic boundary conditions  in z direction
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(false);

	move.init(ingredients);
	move.setPosition(1,-3,-1);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(10,10,0);
	EXPECT_TRUE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,13,15);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,10,10);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,9,9);
	EXPECT_FALSE(move.check(ingredients));

	move.init(ingredients);
	move.setPosition(1,8,8);
	EXPECT_TRUE(move.check(ingredients));




}
/*****************************************************************************/
/**
 * @class FeatureBoxReadTest
 * @brief test fixture for the testing of the reading functions !box_ and !periodic_
 * @details when using TEST_F(FeatureBoxReadTest, -testname-),
 * instances of istringstreams are ready to use \n
 * these structures are used in all the following tests \n
 * the class provides different streams to test the reaction of these
 * eading funcions to different input \n
 * */
/*************************************************************************** */
class FeatureBoxReadTest:public ::testing::Test{
public:
  istringstream inputBoxStandard;
  istringstream inputBoxSpaces;
  istringstream inputBoxMissing;
  istringstream inputBoxText;

  istringstream inputPeriodicStandard;
  istringstream inputPeriodicSpaces;
  istringstream inputPeriodicMissing;
  istringstream inputPeriodicNumber;
  istringstream inputPeriodicText;

  virtual void SetUp(){
    inputBoxStandard.str("100 ");
    inputBoxSpaces.str("    200  ");
    inputBoxMissing.str(" ");
    inputBoxText.str("nonesense");

    inputPeriodicStandard.str("1");
    inputPeriodicSpaces.str("   0 ");
    inputPeriodicMissing.str(" ");
    inputPeriodicNumber.str("4");
    inputPeriodicText.str("nonesense");
  }
};

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxReadTest,ReadBoxX)
 * @brief test the class ReadBoxX
/*************************************************************************** */
TEST_F(FeatureBoxReadTest,ReadBoxX){

  //create ReadBoxX object
  FeatureBox box;
  ReadBoxX<FeatureBox> Read(box);

  //test reaction of the Read object to different streams
  Read.setInputStream(&inputBoxStandard);
  Read.execute();
  EXPECT_EQ(box.getBoxX(),100);

  Read.setInputStream(&inputBoxSpaces);
  Read.execute();
  EXPECT_EQ(box.getBoxX(),200);

  Read.setInputStream(&inputBoxMissing);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputBoxText);
  EXPECT_THROW(Read.execute(),std::runtime_error);

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxReadTest,ReadBoxY)
 * @brief test the class ReadBoxX
/*************************************************************************** */
TEST_F(FeatureBoxReadTest,ReadBoxY){

  //create ReadBoxY object
  FeatureBox box;
  ReadBoxY<FeatureBox> Read(box);

  //test reaction of the Read object to different streams
  Read.setInputStream(&inputBoxStandard);
  Read.execute();
  EXPECT_EQ(box.getBoxY(),100);

  Read.setInputStream(&inputBoxSpaces);
  Read.execute();
  EXPECT_EQ(box.getBoxY(),200);

  Read.setInputStream(&inputBoxMissing);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputBoxText);
  EXPECT_THROW(Read.execute(),std::runtime_error);

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxReadTest,ReadBoxZ)
 * @brief test the class ReadBoxZ
/*************************************************************************** */
TEST_F(FeatureBoxReadTest,ReadBoxZ){

  //create ReadBoxZ object
  FeatureBox box;
  ReadBoxZ<FeatureBox> Read(box);

  //test reaction of the Read object to different streams
  Read.setInputStream(&inputBoxStandard);
  Read.execute();
  EXPECT_EQ(box.getBoxZ(),100);

  Read.setInputStream(&inputBoxSpaces);
  Read.execute();
  EXPECT_EQ(box.getBoxZ(),200);

  Read.setInputStream(&inputBoxMissing);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputBoxText);
  EXPECT_THROW(Read.execute(),std::runtime_error);

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxReadTest,ReadPeriodicX)
 * @brief test the class ReadPeriodicX
/*************************************************************************** */
TEST_F(FeatureBoxReadTest,ReadPeriodicX){

  //create ReadPeriodicX object
  FeatureBox box;
  ReadPeriodicX<FeatureBox> Read(box);

  //test reaction of the Read object to different streams
  Read.setInputStream(&inputPeriodicStandard);
  Read.execute();
  EXPECT_TRUE(box.isPeriodicX());

  Read.setInputStream(&inputPeriodicSpaces);
  Read.execute();
  EXPECT_FALSE(box.isPeriodicX());

  Read.setInputStream(&inputPeriodicMissing);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputPeriodicNumber);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputPeriodicText);
  EXPECT_THROW(Read.execute(),std::runtime_error);

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxReadTest,ReadPeriodicY)
 * @brief test the class ReadPeriodicY
/*************************************************************************** */
TEST_F(FeatureBoxReadTest,ReadPeriodicY){

  //create ReadPeriodicY object
  FeatureBox box;
  ReadPeriodicY<FeatureBox> Read(box);

  //test reaction of the Read object to different streams
  Read.setInputStream(&inputPeriodicStandard);
  Read.execute();
  EXPECT_TRUE(box.isPeriodicY());

  Read.setInputStream(&inputPeriodicSpaces);
  Read.execute();
  EXPECT_FALSE(box.isPeriodicY());

  Read.setInputStream(&inputPeriodicMissing);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputPeriodicNumber);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputPeriodicText);
  EXPECT_THROW(Read.execute(),std::runtime_error);

}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxReadTest,ReadPeriodicZ)
 * @brief test the class ReadPeriodicZ
/*************************************************************************** */
TEST_F(FeatureBoxReadTest,ReadPeriodicZ){

  //create ReadPeriodicZ object
  FeatureBox box;
  ReadPeriodicZ<FeatureBox> Read(box);

  //test reaction of the Read object to different streams
  Read.setInputStream(&inputPeriodicStandard);
  Read.execute();
  EXPECT_TRUE(box.isPeriodicZ());

  Read.setInputStream(&inputPeriodicSpaces);
  Read.execute();
  EXPECT_FALSE(box.isPeriodicZ());

  Read.setInputStream(&inputPeriodicMissing);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputPeriodicNumber);
  EXPECT_THROW(Read.execute(),std::runtime_error);

  Read.setInputStream(&inputPeriodicText);
  EXPECT_THROW(Read.execute(),std::runtime_error);

}

/****************************************************************************/
//test classes FeatureBoxWrite.h
/****************************************************************************/
/*****************************************************************************/
/**
 * @class FeatureBoxTest
 * @brief prepare ingredients system
 * @details when using TEST_F(FeatureBoxTest, -testname-), an instance of
 * Ingredients (ingredients) with features Box and Bondset is prepared. \n
 * for large/less output just comment/uncomment the depending parts of the test classes
 * */
/*************************************************************************** */
class FeatureBoxWriteTest: public ::testing::Test{
protected:
  //stringstream for writing in instead of a file
  ostringstream stream;
  //setup the box
  FeatureBox box;

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
};

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxWriteTest,WriteBoxX)
 * @brief test the class WriteBoxX
/*************************************************************************** */
TEST_F(FeatureBoxWriteTest,WriteBoxX){
  box.setBoxX(10);
  //this is the expected output
  string expected("!box_x=10\n\n");
  //check if all goes well
  WriteBoxX<FeatureBox> writeBox(box);
  writeBox.writeStream(stream);
  string output=stream.str();
  EXPECT_EQ(expected,output);
}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxWriteTest,WriteBoxY)
 * @brief test the class WriteBoxY
/*************************************************************************** */
TEST_F(FeatureBoxWriteTest,WriteBoxY){
  box.setBoxY(10);
  //this is the expected output
  string expected("!box_y=10\n\n");
  //check if all goes well
  WriteBoxY<FeatureBox> writeBox(box);
  writeBox.writeStream(stream);
  string output=stream.str();
  EXPECT_EQ(expected,output);
}

/*****************************************************************************/
/**
 * @fn TEST_F(FeatureBoxWriteTest,WriteBoxZ)
 * @brief test the class WriteBoxZ
/*************************************************************************** */
TEST_F(FeatureBoxWriteTest,WriteBoxZ){
  box.setBoxZ(10);
  //this is the expected output
  string expected("!box_z=10\n\n");
  //check if all goes well
  WriteBoxZ<FeatureBox> writeBox(box);
  writeBox.writeStream(stream);
  string output=stream.str();
  EXPECT_EQ(expected,output);
}

TEST_F(FeatureBoxWriteTest,WritePeriodicX){
  box.setPeriodicX(0);
  //this is the expected output
  string expected("!periodic_x=0\n\n");
  //check if all goes well
  WritePeriodicX<FeatureBox> writeBox(box);
  writeBox.writeStream(stream);
  string output=stream.str();
  EXPECT_EQ(expected,output);
}

TEST_F(FeatureBoxWriteTest,WritePeriodicY){
  box.setPeriodicY(0);
  //this is the expected output
  string expected("!periodic_y=0\n\n");
  //check if all goes well
  WritePeriodicY<FeatureBox> writeBox(box);
  writeBox.writeStream(stream);
  string output=stream.str();
  EXPECT_EQ(expected,output);
}

TEST_F(FeatureBoxWriteTest,WritePeriodicZ){
  box.setPeriodicZ(0);
  //this is the expected output
  string expected("!periodic_z=0\n\n");
  //check if all goes well
  WritePeriodicZ<FeatureBox> writeBox(box);
  writeBox.writeStream(stream);
  string output=stream.str();
  EXPECT_EQ(expected,output);
}

