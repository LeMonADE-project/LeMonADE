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
 * @brief Tests for the classes Molecules
 *
 * @author Christoph
 * @date 23.04.2013
 * */
/*****************************************************************************/

#include <stdexcept>
#include <cstdio>

#include "gtest/gtest.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/core/ConfigureSystem.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class MoleculesTest: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_2(FeatureBox, FeatureBondset<>) Features;
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

TEST_F(MoleculesTest, GetterAndSetter){

	Molecules <VectorInt3> molecules;
    EXPECT_EQ(7, molecules.getMaxConnectivity());
	EXPECT_EQ(0, molecules.getAge());

	//this calles the implicitely defined copy-constructor,
	//not the conversion constructor defined in Molecules.h
	Molecules <VectorInt3> molecules2(molecules);
	EXPECT_EQ(0, molecules2.getAge());

	molecules2.setAge(10000);
	EXPECT_EQ(10000, molecules2.getAge());

	molecules2.resize(10);
	EXPECT_EQ(10, molecules2.size() );

	molecules2.connect(0,1);
	molecules2.connect(1,2,1);


	EXPECT_EQ(1, molecules2.getLinkInfo(1,2));
	EXPECT_EQ(1, molecules2.getLinkInfo(2,1));
	EXPECT_THROW(molecules2.getLinkInfo(1,11),std::runtime_error);



	molecules2.setLinkInfo(1,2, 4);
	EXPECT_EQ(4, molecules2.getLinkInfo(2,1));
	EXPECT_EQ(4, molecules2.getLinkInfo(1,2));

	EXPECT_THROW(molecules2.setLinkInfo(1,11,4),std::runtime_error);

	EXPECT_EQ(2, molecules2.getNumLinks(1));
	EXPECT_EQ(1, molecules2.getNumLinks(2));
	EXPECT_EQ(1, molecules2.getNumLinks(0));
	EXPECT_EQ(0, molecules2.getNumLinks(9));

	EXPECT_THROW(molecules2.getNumLinks(11),std::out_of_range);

	EXPECT_EQ(0, molecules2.getNeighborIdx(1, 0));
	EXPECT_EQ(2, molecules2.getNeighborIdx(1, 1));

}

TEST_F(MoleculesTest, AddMonomers){

	//this function tests the AddMonomer-functionality that returns the monomer idx in the list
	Molecules <VectorInt3> molecules;
    EXPECT_EQ(7, molecules.getMaxConnectivity());
	EXPECT_EQ(0, molecules.getAge());


	EXPECT_EQ(0, molecules.size() );

	EXPECT_EQ(0, molecules.addMonomer(0,4,6));
	EXPECT_EQ(1, molecules.addMonomer(0,6,6));

	molecules.resize(6);
	EXPECT_EQ(6, molecules.size() );
	EXPECT_EQ(6, molecules.addMonomer(0,6,8) );
	EXPECT_EQ(7, molecules.addMonomer(0,6,10) );

	VectorInt3 mono1(2,6,10);
	VectorInt3 mono2(4,6,10);

	EXPECT_EQ(8, molecules.addMonomer(mono1) );
	EXPECT_EQ(9, molecules.addMonomer(mono2) );

	EXPECT_EQ(10, molecules.addMonomer(molecules[9].getX()+2, molecules[9].getY(), molecules[9].getZ() ) );


	//check all coordiantes
	EXPECT_EQ(0, molecules[0].getX());
	EXPECT_EQ(4, molecules[0].getY());
	EXPECT_EQ(6, molecules[0].getZ());

	EXPECT_EQ(0, molecules[1].getX());
	EXPECT_EQ(6, molecules[1].getY());
	EXPECT_EQ(6, molecules[1].getZ());

	for(uint i = 2; i < 6; i++){
		EXPECT_EQ(0, molecules[i].getX());
		EXPECT_EQ(0, molecules[i].getY());
		EXPECT_EQ(0, molecules[i].getZ());
	}

	EXPECT_EQ(0, molecules[6].getX());
	EXPECT_EQ(6, molecules[6].getY());
	EXPECT_EQ(8, molecules[6].getZ());

	EXPECT_EQ(0, molecules[7].getX());
	EXPECT_EQ(6, molecules[7].getY());
	EXPECT_EQ(10, molecules[7].getZ());

	EXPECT_EQ(2, molecules[8].getX());
	EXPECT_EQ(6, molecules[8].getY());
	EXPECT_EQ(10, molecules[8].getZ());

	EXPECT_EQ(4, molecules[9].getX());
	EXPECT_EQ(6, molecules[9].getY());
	EXPECT_EQ(10, molecules[9].getZ());

	EXPECT_EQ(6, molecules[10].getX());
	EXPECT_EQ(6, molecules[10].getY());
	EXPECT_EQ(10, molecules[10].getZ());

	molecules.clear();
	EXPECT_EQ(0, molecules.size() );

}

TEST_F(MoleculesTest, Connectivity){

  //create test object
  Molecules <VectorInt3,3> molecules;
  molecules.resize(10);

  //connect to non-existent monomer
  EXPECT_THROW(molecules.connect(1,11),std::range_error);
  EXPECT_EQ(0,molecules.getNumLinks(1));

  //try to exceed maximum connectivity
  molecules.connect(1,2);
  molecules.connect(1,3);
  molecules.connect(1,4);
  molecules.connect(2,3);
  EXPECT_THROW(molecules.connect(1,5),std::runtime_error);
  EXPECT_EQ(3,molecules.getNumLinks(1));
  EXPECT_EQ(0,molecules.getNumLinks(5));
  EXPECT_EQ(0,molecules.getLinkInfo(1,2));
  EXPECT_EQ(2,molecules.getNeighborIdx(1,0));
  EXPECT_EQ(3,molecules.getNeighborIdx(1,1));
  EXPECT_EQ(4,molecules.getNeighborIdx(1,2));
  EXPECT_ANY_THROW(molecules.getNeighborIdx(1,3));
  EXPECT_TRUE(molecules.areConnected(1,2));
  EXPECT_FALSE(molecules.areConnected(2,4));
  EXPECT_FALSE(molecules.areConnected(1,20));

  //disconnect monomers
  molecules.disconnect(1,2);
  EXPECT_EQ(2,molecules.getNumLinks(1));
  EXPECT_EQ(1,molecules.getNumLinks(2));
  EXPECT_ANY_THROW(molecules.getLinkInfo(1,2));

  //try to disconnect monomers which are not connected
  EXPECT_ANY_THROW(molecules.disconnect(1,2));
}

TEST_F(MoleculesTest, Constructors){

  //standard constructor
  Molecules <VectorInt3> molecules1;
  //check default values
  EXPECT_EQ(0,molecules1.getAge());
  EXPECT_EQ(0,molecules1.size());

  //conversion constructor
  //add some information to molecules 1 and see if the conversion constructor
  //does a good copying job
  molecules1.resize(10);
  molecules1.connect(0,1);
  molecules1.connect(0,2);
  molecules1.connect(0,3);
  molecules1.connect(0,4,2);
  molecules1[1].setAllCoordinates(2,3,4);
  //this should throw an exception, because the max connectivity of molecules1
  //is too high
  EXPECT_ANY_THROW( (Molecules<VectorInt3, 2> (molecules1)) );

  //create new objects from molecules1
  Molecules<VectorInt3> molecules2(molecules1);
  Molecules<VectorInt3, 8> molecules3(molecules2);

  //check if everything was copied correctly
  EXPECT_EQ(7,molecules2.getMaxConnectivity());
  EXPECT_EQ(8,molecules3.getMaxConnectivity());
  EXPECT_EQ(10,molecules2.size());
  EXPECT_EQ(10,molecules3.size());
  EXPECT_EQ(molecules1[1],molecules2[1]);
  EXPECT_EQ(molecules2[1],molecules3[1]);
  EXPECT_EQ(molecules1.getLinkInfo(0,4),molecules2.getLinkInfo(0,4));
  EXPECT_EQ(molecules2.getLinkInfo(0,4),molecules3.getLinkInfo(0,4));
}

TEST_F(MoleculesTest, AssignmentOperator){

  //setup of a molecules object first
  Molecules <VectorInt3> molecules1;

  molecules1.resize(10);
  molecules1.connect(0,1);
  molecules1.connect(0,2);
  molecules1.connect(0,3);
  molecules1.connect(0,4,2);
  molecules1[1].setAllCoordinates(2,3,4);

  //copy to molecules object of same type, and with higher max
  //connectivity...should both work
  Molecules<VectorInt3> molecules2;
  molecules2=molecules1;
  Molecules<VectorInt3, 8> molecules3;
  molecules3=molecules2;

  //check if everything was copied correctly
  EXPECT_EQ(7,molecules2.getMaxConnectivity());
  EXPECT_EQ(8,molecules3.getMaxConnectivity());
  EXPECT_EQ(10,molecules2.size());
  EXPECT_EQ(10,molecules3.size());
  EXPECT_EQ(molecules1[1],molecules2[1]);
  EXPECT_EQ(molecules2[1],molecules3[1]);
  EXPECT_EQ(molecules1.getLinkInfo(0,4),molecules2.getLinkInfo(0,4));
  EXPECT_EQ(molecules2.getLinkInfo(0,4),molecules3.getLinkInfo(0,4));


  //now try to assign to molecules object with too low max connectivity
  //this should throw an exception
  Molecules<VectorInt3, 2> molecules4;
  EXPECT_ANY_THROW( molecules4=molecules1) ;

}

TEST_F(MoleculesTest,Clear){
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

  EXPECT_EQ(8,ingredients.getMolecules().getTotalNumLinks());
  ingredients.modifyMolecules().clearBonds();
  EXPECT_EQ(0,ingredients.getMolecules().getTotalNumLinks());
  for(size_t n=0;n<ingredients.getMolecules().size();n++)
  {
	  EXPECT_EQ(0,ingredients.getMolecules().getNumLinks(n));
  }

  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().connect(3,4);
  ingredients.modifyMolecules().connect(4,5);
  ingredients.modifyMolecules().connect(6,7);
  ingredients.modifyMolecules().connect(7,8);
  ingredients.modifyMolecules().connect(8,9);
  ingredients.modifyMolecules().connect(1,4);

  EXPECT_EQ(8,ingredients.getMolecules().getTotalNumLinks());
  ingredients.modifyMolecules().clear();
  EXPECT_EQ(0,ingredients.getMolecules().getTotalNumLinks());
  EXPECT_EQ(0,ingredients.getMolecules().size());

}
