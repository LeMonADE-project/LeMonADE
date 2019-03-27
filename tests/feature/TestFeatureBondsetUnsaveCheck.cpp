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
#include <LeMonADE/feature/FeatureBondsetUnsaveCheck.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/GenerateMonomerType.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveConnectSc.h>

using namespace std;
/************************************************************************/
//define test fixtures for the different tests their purpose is to set up
//the tests to suppress cout's output such that is does not display on the
//standard output during the tests. this makes google test's output more readeable
/************************************************************************/

class FeatureBondsetTestUnsaveCheck: public ::testing::Test{
public:
  typedef LOKI_TYPELIST_3(FeatureConnectionSc, FeatureBondsetUnsaveCheck<>,FeatureBox) Features;
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

/***********************************************************************/
//checkMove for unknown moves should return true
/***********************************************************************/
TEST_F(FeatureBondsetTestUnsaveCheck,checkUnknownMove)
{
	UnknownMove someMove;
	FeatureBondset<> feature;
	EXPECT_TRUE(someMove.check(feature));
}

/***********************************************************************/
//checkMove for MoveLocalBase
//same check as for FeatureBondset.
//The FeatureBondsetUnsaveCheck should redirect the checkMove() to the 
//FeatureBondset
/***********************************************************************/
TEST_F(FeatureBondsetTestUnsaveCheck,checkLocalScBfmMove)
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

/***********************************************************************/
//checkMove for MoveLocalBase
/***********************************************************************/
TEST_F(FeatureBondsetTestUnsaveCheck,checkConnectmMoveSc)
{
	Ing ingredients;

	ingredients.setPeriodicX(true);
	ingredients.setPeriodicY(true);
	ingredients.setPeriodicZ(true);
	ingredients.setBoxX(16);
	ingredients.setBoxY(16);
	ingredients.setBoxZ(16);

	ingredients.modifyMolecules().resize(2);
	ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
	ingredients.modifyMolecules()[1].setAllCoordinates(18,0,0);

	ingredients.modifyMolecules()[0].setReactive(true);
	ingredients.modifyMolecules()[1].setReactive(true);
	ingredients.modifyMolecules()[0].setNumMaxLinks(1);
	ingredients.modifyMolecules()[1].setNumMaxLinks(1);
	EXPECT_TRUE (ingredients.getMolecules()[0].isReactive() );
	EXPECT_TRUE(ingredients.getMolecules()[1].isReactive() );
	EXPECT_EQ(ingredients.getMolecules()[0].getNumMaxLinks(),1);
	EXPECT_EQ(ingredients.getMolecules()[1].getNumMaxLinks(),1);

	ingredients.modifyBondset().addBond(2,0,0,17);
	ingredients.modifyBondset().addBond(-2,0,0,20);

	ingredients.synchronize(ingredients);
	EXPECT_EQ(0,ingredients.getIdFromLattice(0,0,0));    
	EXPECT_EQ(1,ingredients.getIdFromLattice(2,0,0));    
	MoveConnectSc move; 
	move.init(ingredients,0,VectorInt3(2,0,0));
	EXPECT_EQ(1,move.getPartner());    
	EXPECT_EQ(0,move.getIndex());    
	EXPECT_EQ(VectorInt3(2,0,0),move.getDir());   
	VectorInt3 bond(ingredients.getMolecules()[1]-ingredients.getMolecules()[0]);
	EXPECT_EQ(VectorInt3(18,0,0),bond);
	//refold bond vector with &7 (only applicable for box dimensions of power 2)
	bond.setAllCoordinates(bond.getX() & 7, bond.getY() & 7, bond.getZ() & 7 );
	EXPECT_EQ(VectorInt3(2,0,0),bond);
	EXPECT_TRUE( ingredients.getBondset().isValid(bond) );
	EXPECT_TRUE(move.check(ingredients));
	
}


