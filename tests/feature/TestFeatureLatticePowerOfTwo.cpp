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

#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureLatticePowerOfTwo.h>

using namespace std;
/*****************************************************************************/
/**
 * @file
 * @brief Tests for the FeatureLatticePowerOfTwo
 *
 * @author Ron
 * @date 17.07.2014
 * */
/*****************************************************************************/

/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class FeatureLatticePowerOfTwoTest: public ::testing::Test{
public:
  //typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
 // typedef ConfigureSystem<VectorInt3,Features> Config;
  //typedef Ingredients<Config> Ing;

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
//checks if synchronize works and throws exceptions for lattice type bool
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, Synchronize1){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(7);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);


	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(129);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(128);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty
	for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
					EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (bool) 0);


}
/************************************************************************/
//checks if synchronize works and throws exceptions for lattice type uint8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, Synchronize2){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(7);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);


	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(129);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(128);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty
	for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
					EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (uint8_t) 0);


}


/************************************************************************/
//checks if synchronize works and throws exceptions for lattice type int8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, Synchronize3){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint32_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(7);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);


	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(129);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(128);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty
	for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
					EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (int8_t) 0);


}

/************************************************************************/
//checks if synchronize works and throws exceptions for lattice type uint32_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, Synchronize4){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint32_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(7);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);


	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(129);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(128);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty
	for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
					EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (uint32_t) 0);


}


/************************************************************************/
//checks if synchronize works and throws exceptions for lattice type int32_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, Synchronize5){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint32_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(7);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);


	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(255);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(129);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(129);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_THROW(ingredients.synchronize(ingredients),std::runtime_error);

	ingredients.setBoxX(128);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(8);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty
	for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
					EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (int32_t) 0);


}


/************************************************************************/
//checks if set and getLatticeEntry for lattice type bool
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice1){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (bool) 0);

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(x,y,z, true);

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

}


/************************************************************************/
//checks if set and getLatticeEntry for lattice type uint8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice2){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), uint8_t (0));

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(x,y,z, (uint8_t) 1);

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(x,y,z, (uint8_t) -1);

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

		EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*((uint8_t) -1));

}

/************************************************************************/
//checks if set and getLatticeEntry for lattice type int8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice3){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<int8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), int8_t (0));

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(x,y,z, (int8_t) 1);

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(x,y,z, (int8_t) -1);

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

		EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*((int8_t) -1));

}

/************************************************************************/
//checks if set and getLatticeEntry for lattice type bool
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice4){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (bool) 0);

	VectorInt3 id(0,0,0);

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
			{
				id.setAllCoordinates(x,y,z);
				ingredients.setLatticeEntry(id, true);

			}

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

}


/************************************************************************/
//checks if set and getLatticeEntry for lattice type uint8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice5){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), uint8_t (0));

	VectorInt3 id(0,0,0);
	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
			{
				id.setAllCoordinates(x,y,z);
				ingredients.setLatticeEntry(id, (uint8_t) 1);
			}

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
			{
				id.setAllCoordinates(x,y,z);
				ingredients.setLatticeEntry(id, (uint8_t) -1);
			}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

		EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*((uint8_t) -1));

}

/************************************************************************/
//checks if set and getLatticeEntry for lattice type int8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice6){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<int8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), int8_t (0));

	VectorInt3 id(0,0,0);
	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
			{
				id.setAllCoordinates(x,y,z);
				ingredients.setLatticeEntry(id, (int8_t) 1);
			}

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
			{
				id.setAllCoordinates(x,y,z);
				ingredients.setLatticeEntry(id, (int8_t) -1);
			}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

		EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*((int8_t) -1));

}



/************************************************************************/
//checks if set and getLatticeEntry and folding for lattice type bool
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice7){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (bool) 0);

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(-5006375+x,42324+y,-64263+z, true);

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());


	// clear lattice
		for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
					ingredients.setLatticeEntry(x,y,z, false);



	VectorInt3 vi(-5006375,42324,-64263);

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(vi+VectorInt3(x,y,z), true);

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());

}


/************************************************************************/
//checks if set and getLatticeEntry and folding for lattice type uint8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice8){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), uint8_t (0));

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(-5006375+x,42324+y,-64263+z, uint8_t (-1));

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*(uint8_t (-1)));

	// clear lattice
		for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
						ingredients.setLatticeEntry(x,y,z, false);


	VectorInt3 vi(-5006375,42324,-64263);

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(vi+VectorInt3(x,y,z), uint8_t (-2));

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*(uint8_t (-2)));

}

/************************************************************************/
//checks if set and getLatticeEntry and folding for lattice type int8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, SetGetLattice9){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<int8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), int8_t (0));

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(-5006375+x,42324+y,-64263+z, int8_t (-1));

	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*(int8_t (-1)));

	// clear lattice
		for(int x=0; x < ingredients.getBoxX(); x++)
			for(int y=0; y < ingredients.getBoxY(); y++)
				for(int z=0; z < ingredients.getBoxZ(); z++)
						ingredients.setLatticeEntry(x,y,z, false);


	VectorInt3 vi(-5006375,42324,-64263);

	// fill all places with true
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				ingredients.setLatticeEntry(vi+VectorInt3(x,y,z), int8_t (-2));

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ()*(int8_t (-2)));

}









/************************************************************************/
//checks the moveOnLattice for lattice type bool
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, MoveOnLattice1){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (bool) 0);

	// fill all places with true
	ingredients.setLatticeEntry(0,0,0, true);

	VectorInt3 viOld(0,0,0);
	VectorInt3 viNew(0,0,0);

	// move first oldPlace to Newlace
	for(int x=1; x < ingredients.getBoxX(); x++)
		{
			viNew += VectorInt3(1,0,0);
			ingredients.moveOnLattice(viOld, viNew);
			viOld += VectorInt3(1,0,0);
		}


	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, 1);
	EXPECT_EQ(ingredients.getLatticeEntry(0,0,0), false);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), true);

	// move first oldPlace to Newlace
	for(int y=1; y < ingredients.getBoxY(); y++)
	{
		viNew += VectorInt3(0,1,0);
		ingredients.moveOnLattice(viOld, viNew);
		viOld += VectorInt3(0,1,0);
	}


	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, 1);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), false);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), true);


	// move first oldPlace to Newlace
	for(int z=1; z < ingredients.getBoxZ(); z++)
		{
			viNew += VectorInt3(0,0,1);
			ingredients.moveOnLattice(viOld, viNew);
			viOld += VectorInt3(0,0,1);
		}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, 1);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), false);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,ingredients.getBoxZ()-1), true);

}

/************************************************************************/
//checks the moveOnLattice for lattice type uint8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, MoveOnLattice2){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (uint8_t) 0);

	// fill all places with true
	ingredients.setLatticeEntry(0,0,0, uint8_t (-1));

	VectorInt3 viOld(0,0,0);
	VectorInt3 viNew(0,0,0);

	// move first oldPlace to Newlace
	for(int x=1; x < ingredients.getBoxX(); x++)
		{
			viNew += VectorInt3(1,0,0);
			ingredients.moveOnLattice(viOld, viNew);
			viOld += VectorInt3(1,0,0);
		}


	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, uint8_t (-1));
	EXPECT_EQ(ingredients.getLatticeEntry(0,0,0), uint8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), uint8_t (-1));

	// move first oldPlace to Newlace
	for(int y=1; y < ingredients.getBoxY(); y++)
	{
		viNew += VectorInt3(0,1,0);
		ingredients.moveOnLattice(viOld, viNew);
		viOld += VectorInt3(0,1,0);
	}


	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, uint8_t (-1));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), uint8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), uint8_t (-1));


	// move first oldPlace to Newlace
	for(int z=1; z < ingredients.getBoxZ(); z++)
		{
			viNew += VectorInt3(0,0,1);
			ingredients.moveOnLattice(viOld, viNew);
			viOld += VectorInt3(0,0,1);
		}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, uint8_t (-1));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), uint8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,ingredients.getBoxZ()-1), uint8_t (-1));

}

/************************************************************************/
//checks the moveOnLattice for lattice type int8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, MoveOnLattice3){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<int8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (int8_t) 0);

	// fill all places with true
	ingredients.setLatticeEntry(0,0,0, int8_t (-127));

	VectorInt3 viOld(0,0,0);
	VectorInt3 viNew(0,0,0);

	// move first oldPlace to Newlace
	for(int x=1; x < ingredients.getBoxX(); x++)
		{
			viNew += VectorInt3(1,0,0);
			ingredients.moveOnLattice(viOld, viNew);
			viOld += VectorInt3(1,0,0);
		}


	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, int8_t (-127));
	EXPECT_EQ(ingredients.getLatticeEntry(0,0,0), int8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), int8_t (-127));

	// move first oldPlace to Newlace
	for(int y=1; y < ingredients.getBoxY(); y++)
	{
		viNew += VectorInt3(0,1,0);
		ingredients.moveOnLattice(viOld, viNew);
		viOld += VectorInt3(0,1,0);
	}


	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, int8_t (-127));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), int8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), int8_t (-127));


	// move first oldPlace to Newlace
	for(int z=1; z < ingredients.getBoxZ(); z++)
		{
			viNew += VectorInt3(0,0,1);
			ingredients.moveOnLattice(viOld, viNew);
			viOld += VectorInt3(0,0,1);
		}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, int8_t (-127));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), int8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,ingredients.getBoxZ()-1), int8_t (-127));

}


/************************************************************************/
//checks the moveOnLattice for lattice type bool
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, MoveOnLattice4){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<bool> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(128);
	ingredients.setBoxY(32);
	ingredients.setBoxZ(64);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (bool) 0);

	// fill all places with true
	ingredients.setLatticeEntry(0,0,0, true);

	int viOld= 0;
	int viNew= 0;

	// move first oldPlace to Newlace
	for(int x=1; x < ingredients.getBoxX(); x++)
		{
			viNew += 1;
			ingredients.moveOnLattice(viOld, 0, 0, viNew,0, 0);
			viOld += 1;
		}


	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, 1);
	EXPECT_EQ(ingredients.getLatticeEntry(0,0,0), false);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), true);

	viOld= 0;
	viNew= 0;

	// move first oldPlace to Newlace
	for(int y=1; y < ingredients.getBoxY(); y++)
	{
		viNew += 1;
		ingredients.moveOnLattice(ingredients.getBoxX()-1,viOld,0, ingredients.getBoxX()-1, viNew, 0);
		viOld += 1;
	}


	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, 1);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), false);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), true);


	viOld= 0;
	viNew= 0;
	// move first oldPlace to Newlace
	for(int z=1; z < ingredients.getBoxZ(); z++)
		{
			viNew += 1;
			ingredients.moveOnLattice(ingredients.getBoxX()-1,ingredients.getBoxY()-1, viOld, ingredients.getBoxX()-1,ingredients.getBoxY()-1, viNew);
			viOld += 1;
		}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z)?1:0;

	EXPECT_EQ(counter, 1);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), false);
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,ingredients.getBoxZ()-1), true);

}

/************************************************************************/
//checks the moveOnLattice for lattice type uint8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, MoveOnLattice5){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<uint8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (uint8_t) 0);

	// fill all places with true
	ingredients.setLatticeEntry(0,0,0, uint8_t (-1));

	int viOld= 0;
	int viNew= 0;

	// move first oldPlace to Newlace
	for(int x=1; x < ingredients.getBoxX(); x++)
		{
			viNew += 1;
			ingredients.moveOnLattice(viOld, 0, 0, viNew,0, 0);
			viOld += 1;
		}


	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, uint8_t (-1));
	EXPECT_EQ(ingredients.getLatticeEntry(0,0,0), uint8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), uint8_t (-1));

	viOld= 0;
	viNew= 0;

	// move first oldPlace to Newlace
	for(int y=1; y < ingredients.getBoxY(); y++)
	{
		viNew += 1;
		ingredients.moveOnLattice(ingredients.getBoxX()-1,viOld,0, ingredients.getBoxX()-1, viNew, 0);
		viOld += 1;
	}


	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, uint8_t (-1));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), uint8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), uint8_t (-1));


	viOld= 0;
	viNew= 0;
	// move first oldPlace to Newlace
	for(int z=1; z < ingredients.getBoxZ(); z++)
		{
			viNew += 1;
			ingredients.moveOnLattice(ingredients.getBoxX()-1,ingredients.getBoxY()-1, viOld, ingredients.getBoxX()-1,ingredients.getBoxY()-1, viNew);
			viOld += 1;
		}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, uint8_t (-1));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), uint8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,ingredients.getBoxZ()-1), uint8_t (-1));

}


/************************************************************************/
//checks the moveOnLattice for lattice type int8_t
/************************************************************************/
TEST_F(FeatureLatticePowerOfTwoTest, MoveOnLattice6){

	typedef LOKI_TYPELIST_1(FeatureLatticePowerOfTwo<int8_t> ) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;

	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(8);
	ingredients.setBoxY(256);
	ingredients.setBoxZ(512);

	EXPECT_NO_THROW(ingredients.synchronize(ingredients));

	// check if all places are empty (false)
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				EXPECT_EQ(ingredients.getLatticeEntry(x,y,z), (int8_t) 0);

	// fill all places with true
	ingredients.setLatticeEntry(0,0,0, int8_t (-255));

	int viOld= 0;
	int viNew= 0;

	// move first oldPlace to Newlace
	for(int x=1; x < ingredients.getBoxX(); x++)
		{
			viNew += 1;
			ingredients.moveOnLattice(viOld, 0, 0, viNew,0, 0);
			viOld += 1;
		}


	//count the lattice entries
	uint32_t counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, int8_t (-255));
	EXPECT_EQ(ingredients.getLatticeEntry(0,0,0), int8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), int8_t (-255));

	viOld= 0;
	viNew= 0;

	// move first oldPlace to Newlace
	for(int y=1; y < ingredients.getBoxY(); y++)
	{
		viNew += 1;
		ingredients.moveOnLattice(ingredients.getBoxX()-1,viOld,0, ingredients.getBoxX()-1, viNew, 0);
		viOld += 1;
	}


	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, int8_t (-255));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,0,0), int8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), int8_t (-255));


	viOld= 0;
	viNew= 0;
	// move first oldPlace to Newlace
	for(int z=1; z < ingredients.getBoxZ(); z++)
		{
			viNew += 1;
			ingredients.moveOnLattice(ingredients.getBoxX()-1,ingredients.getBoxY()-1, viOld, ingredients.getBoxX()-1,ingredients.getBoxY()-1, viNew);
			viOld += 1;
		}

	//count the lattice entries
	counter=0;
	for(int x=0; x < ingredients.getBoxX(); x++)
		for(int y=0; y < ingredients.getBoxY(); y++)
			for(int z=0; z < ingredients.getBoxZ(); z++)
				counter+=ingredients.getLatticeEntry(x,y,z);

	EXPECT_EQ(counter, int8_t (-255));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,0), int8_t (0));
	EXPECT_EQ(ingredients.getLatticeEntry(ingredients.getBoxX()-1,ingredients.getBoxY()-1,ingredients.getBoxZ()-1), int8_t (-255));

}
