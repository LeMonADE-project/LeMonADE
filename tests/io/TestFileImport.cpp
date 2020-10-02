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

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/io/AbstractRead.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class FileImportTest: public ::testing::Test{
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

/*
 * this test case only creates an object of the class FileImport
 * and checks if everything that is supposed to be done in the constructor
 * works correctly.
 * */
TEST_F(FileImportTest,Initialization)
{
  //create a system using FeatureBondset and the Molecules template
  typedef LOKI_TYPELIST_1(FeatureMoleculesIO)	Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients < Config> MyIngredients;
  MyIngredients ingredients;


  //try to open a non-existent file, which should result in an exception
  EXPECT_THROW(FileImport<MyIngredients> nonExistentFile("nofile.bfm",ingredients) ,std::runtime_error );

  //Import the file using the class FileImport
  FileImport<MyIngredients> file("tests/fileImportTest.test",ingredients);

  //scan file for !mcs and read-in first frame
  file.initialize();

  //check interface getDestination()
  EXPECT_EQ(&ingredients, &file.getDestination());
  //check if the ingredients got the correct name
  EXPECT_EQ("tests/fileImportTest.test",ingredients.getName());

  //check if header and connections were read correctly on initializing the FileImport

  //check number of monomers
  EXPECT_EQ(10,ingredients.getMolecules().size());

  //check some bond information:
  //number of bonds
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(0));
  EXPECT_EQ(2, ingredients.getMolecules().getNumLinks(3));
  EXPECT_EQ(0, ingredients.getMolecules().getNumLinks(9));
  //bond partners
  EXPECT_EQ(0, ingredients.getMolecules().getNeighborIdx(5,0));
  EXPECT_EQ(5, ingredients.getMolecules().getNeighborIdx(0,0));
  EXPECT_EQ(1, ingredients.getMolecules().getNeighborIdx(0,1));
  //box
  EXPECT_EQ(234,ingredients.getBoxX());
  EXPECT_EQ(453,ingredients.getBoxY());
  EXPECT_EQ(67,ingredients.getBoxZ());
  //periodicity
  EXPECT_TRUE(ingredients.isPeriodicX());
  EXPECT_FALSE(ingredients.isPeriodicY());
  EXPECT_TRUE(ingredients.isPeriodicZ());





}

/* *******************************************
 * in this test the function read() is tested.
 * *******************************************/

TEST_F(FileImportTest, ReadNextMcs)
{
  //create a system using FeatureBondset and the Molecules template
  typedef LOKI_TYPELIST_1(FeatureMoleculesIO)	Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients < Config> MyIngredients;
  MyIngredients ingredients;

  //Import the file using the class FileImport
  FileImport<MyIngredients> file("tests/fileImportTest.test",ingredients);

  //scan file for !mcs and read-in first frame
  file.initialize();

  //read fist mcs. the read function should return true, because the file
  //has not reached the end
  EXPECT_TRUE(file.read());

  //check some positions
  VectorInt3 position2; position2.setAllCoordinates(5,5,5);
  VectorInt3 position3; position3.setAllCoordinates(3,5,4);
  VectorInt3 position4; position4.setAllCoordinates(5,7,5);
  EXPECT_EQ(position2,ingredients.getMolecules()[2]);
  EXPECT_EQ(position3,ingredients.getMolecules()[3]);
  EXPECT_EQ(position4,ingredients.getMolecules()[4]);

  EXPECT_TRUE(file.read());
  EXPECT_TRUE(file.read());
  //check positions and age again
  EXPECT_EQ(30,ingredients.getMolecules().getAge());
  EXPECT_EQ(position2,ingredients.getMolecules()[2]);
  EXPECT_EQ(position3,ingredients.getMolecules()[3]);
  EXPECT_EQ(position4,ingredients.getMolecules()[4]);
  VectorInt3 position9; position9.setAllCoordinates(20,20,20);
  EXPECT_EQ(position9,ingredients.getMolecules()[9]);
  //now there is no mcs left and the function should return false
  EXPECT_FALSE(file.read());

}

/* *******************************************
 * test the function registerRead()
 * for this case, make a dummy command that reads a word
 * from a the in put stream. on the way, the test also
 * checks if the interface close() works.
 * ********************************************/

//dummy command, who's execute() method reads a word from a stream
template<class Destination>
class DummyCommand: public ReadToDestination<Destination>
{
public:
  DummyCommand (Destination& dst):ReadToDestination<Destination>(dst){}
  void execute()
  {
    string aWord;
    this->getInputStream()>>aWord;
    (this->getDestination().contents)+=aWord;
  }
};


//now define three dummy features, that register
//the dummy command.
//the first one registers only one command
class DummyFeature1:public Feature
{
public:
  string contents;
  template<class File> void exportRead(File& file)
  {
    file.registerRead("!firstCommand",new DummyCommand<DummyFeature1>(*this));
  };
  template<class Data> void synchronize(Data&){};
  void setName(string name){};

  // please ignore this function, this was an akward quick fix of a compilation problem
   void addComment(std::string comment){return;}
};

//the second one registers two commands with different names
class DummyFeature2:public Feature
{
public:
  string contents;
  template<class File> void exportRead(File& file)
  {
    file.registerRead("!firstCommand",new DummyCommand<DummyFeature2>(*this));
    file.registerRead("!secondCommand",new DummyCommand<DummyFeature2>(*this));
  };
  template<class Data> void synchronize(Data&){};
  void setName(string name){};

  // please ignore this function, this was an akward quick fix of a compilation problem
   void addComment(std::string comment){return;}
};

//the third one tries to register two commands with the same name,
//which should fail
class DummyFeature3:public Feature
{
public:
  string contents;
  template<class File> void exportRead(File& file)
  {
    file.registerRead("!firstCommand",new DummyCommand<DummyFeature3>(*this));
    file.registerRead("!firstCommand",new DummyCommand<DummyFeature3>(*this));
  };
  template<class Data> void synchronize(Data&){};
  void setName(string name){};

  // please ignore this function, this was an akward quick fix of a compilation problem
   void addComment(std::string comment){return;}
};


//this test case uses the above defined commands and ingredients dummies
//to read the string "freshLemonade" from a file

TEST_F(FileImportTest, RegisterRead)
{
  typedef LOKI_TYPELIST_3(FeatureBox,FeatureBondset<>,DummyFeature1) Features1;
  typedef ConfigureSystem<VectorInt3,Features1> Config1;
  typedef Ingredients<Config1> MyIngredients1;
  MyIngredients1 ingredients1;

  typedef LOKI_TYPELIST_3(FeatureBox,FeatureBondset<>,DummyFeature2) Features2;
  typedef ConfigureSystem<VectorInt3,Features2> Config2;
  typedef Ingredients<Config2> MyIngredients2;
  MyIngredients2 ingredients2;

  typedef LOKI_TYPELIST_3(FeatureBox,FeatureBondset<>,DummyFeature3) Features3;
  typedef ConfigureSystem<VectorInt3,Features3> Config3;
  typedef Ingredients<Config3> MyIngredients3;
  MyIngredients3 ingredients3;


  //Import the file using the class FileImport
  //here we use the first ingredients class, that only knows the command
  // !firstCommand
  FileImport<MyIngredients1> file1("tests/fileImportTest2.test",ingredients1);

  //scan file for !mcs and read-in first frame
  file1.initialize();


  EXPECT_EQ("fresh",ingredients1.contents);
  file1.close();

  //here we use the second ingredients class, which knows both of the commands
  //used in the testfile
  FileImport<MyIngredients2> file2("tests/fileImportTest2.test",ingredients2);

  //scan file for !mcs and read-in first frame
  file2.initialize();


  EXPECT_EQ("freshLemonade",ingredients2.contents);
  file1.close();

  //here we use the third ingredients class, which tries to register two commands
  //with the same command string.
  EXPECT_THROW(FileImport<MyIngredients3> file3("tests/fileImportTest2.test",ingredients3),
	       std::runtime_error);
}

/* *****************************************************************************
 * check the routines ScanFile, gotoMcs, gotoEnd, gotoStart
 */

TEST_F(FileImportTest, gotoMcs)
{
	typedef LOKI_TYPELIST_1(FeatureMoleculesIO)	Features;
	typedef ConfigureSystem<VectorInt3,Features,5> Config;
	typedef Ingredients < Config> MyIngredients;
	MyIngredients ingredients;
	MyIngredients ingredients2;
	//Import the file using the class FileImport
	FileImport<MyIngredients> file("tests/fileImportTest3.test",ingredients);

	//scan file for !mcs and read-in first frame
	file.initialize();

	//save the first conformation (mcs=1000), the last conformation(mcs=2000000)
	//and conformations at mcs=64000 and 236000

	MyIngredients::molecules_type mol1000; 		// first conformation
	MyIngredients::molecules_type mol2000; 		// second conformation
	MyIngredients::molecules_type mol64000;		// 32nd conformation
	MyIngredients::molecules_type mol236000;	// 118th conformation
	MyIngredients::molecules_type mol2000000;	// 1000th conformation



	for(size_t n=0;n<=1000;n++)
	{
		file.read();
		if(n==0) mol1000=ingredients.getMolecules();
		if(n==1) mol2000=ingredients.getMolecules();
		if(n==32) mol64000=ingredients.getMolecules();
		if(n==118) mol236000=ingredients.getMolecules();
		if(n==1000) mol2000000=ingredients.getMolecules();
	}

	file.close();


	//now read the same file into ingredients2 and use the goto routines
	FileImport<MyIngredients> file2("tests/fileImportTest3.test",ingredients2);

	//scan file for !mcs and read-in first frame
	file2.initialize();


	file2.gotoMcs(236000);
	EXPECT_EQ(ingredients2.getMolecules().getAge(),236000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol236000[n]);
	}

 	file2.gotoMcs(1000);
 	EXPECT_EQ(ingredients2.getMolecules().getAge(),1000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol1000[n]);
	}

	file2.read();
	EXPECT_EQ(ingredients2.getMolecules().getAge(),2000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol2000[n]);
	}
//
 	file2.gotoMcs(2000000);
 	EXPECT_EQ(ingredients2.getMolecules().getAge(),2000000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol2000000[n]);
	}

 	file2.gotoMcs(64000);
 	EXPECT_EQ(ingredients2.getMolecules().getAge(),64000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol64000[n]);
	}

 	file2.gotoMcs(236001);
 	EXPECT_EQ(ingredients2.getMolecules().getAge(),238000);

 	file2.gotoMcs(0);
 	EXPECT_EQ(ingredients2.getMolecules().getAge(),1000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol1000[n]);
	}

 	file2.gotoMcs(3000000);
 	EXPECT_EQ(ingredients2.getMolecules().getAge(),2000000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol2000000[n]);
	}

	file2.gotoStart();
	EXPECT_EQ(ingredients2.getMolecules().getAge(),1000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol1000[n]);
	}

	file2.gotoEnd();
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol2000000[n]);
	}

	file2.gotoStart();
	EXPECT_EQ(ingredients2.getMolecules().getAge(),1000);
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol1000[n]);
	}
}


/* *****************************************************************************
 * check the routines ScanFile, gotoMcsSave, gotoEndSave
 */

TEST_F(FileImportTest, gotoMcsSave)
{
	typedef LOKI_TYPELIST_1(FeatureMoleculesIO)	Features;
	typedef ConfigureSystem<VectorInt3,Features,5> Config;
	typedef Ingredients < Config> MyIngredients;
	MyIngredients ingredients;
	MyIngredients ingredients2;
	//Import the file using the class FileImport
	FileImport<MyIngredients> file("tests/fileImportTest3.test",ingredients);

	//scan file for !mcs and read-in first frame
	file.initialize();



	//save the last configuration in the file. go there step by step
	MyIngredients::molecules_type mol2000000;	// 1000th conformation

	for(size_t n=0;n<=1000;n++)
	{
		file.read();
		if(n==1000) mol2000000=ingredients.getMolecules();
	}

	file.close();


	//now read the same file into ingredients2 and use the gotoEndSave
	FileImport<MyIngredients> file2("tests/fileImportTest3.test",ingredients2);

	//scan file for !mcs and read-in first frame
	file2.initialize();


	//check if all positions agree
	file2.gotoEndSave();
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{
		EXPECT_EQ(ingredients2.getMolecules()[n],mol2000000[n]);
	}

	//also check all bonds in this case
	for(size_t n=0;n<ingredients2.getMolecules().size();n++)
	{

		EXPECT_EQ(ingredients2.getMolecules().getNumLinks(n),mol2000000.getNumLinks(n));
		for(size_t m=0;m<n;m++)
		{
			EXPECT_EQ(ingredients2.getMolecules().areConnected(m,n),mol2000000.areConnected(m,n));
		}
	}

}

/* *****************************************************************************
 * check the routines getMinAge(), getMaxAge()
 */

TEST_F(FileImportTest, getMinMaxAge)
{
    typedef LOKI_TYPELIST_1(FeatureMoleculesIO)	Features;
    typedef ConfigureSystem<VectorInt3,Features,5> Config;
    typedef Ingredients < Config> MyIngredients;
    MyIngredients ingredients;
    //Import the file using the class FileImport
    FileImport<MyIngredients> file("tests/fileImportTest3.test",ingredients);

    //scan file for !mcs and read-in first frame
    file.initialize();


    EXPECT_EQ(2000000,file.getMaxAge());
    EXPECT_EQ(1000,file.getMinAge());
}
