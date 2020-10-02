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

/*
 * TestIngedients.cpp
 *
 *  Created on: 31.07.2013
 *      Author: christoph
 */

#include "gtest/gtest.h"

#include <stdexcept>
#include <cstdio>
#include <sstream>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/core/ConfigureSystem.h>

using namespace std;

TEST(IngredientsTest, Constructors){

	typedef LOKI_TYPELIST_2(FeatureBox, FeatureBondset<>) Features;
	typedef ConfigureSystem<VectorInt3,Features,3> Config;
	typedef Ingredients < Config> MyIngredients;
	MyIngredients ingredientsWithName("myName");

	EXPECT_STREQ("myName",ingredientsWithName.getName().c_str());

	MyIngredients ingredientsWithOutName;
	EXPECT_STREQ("new_lemonade" ,ingredientsWithOutName.getName().c_str());
}

TEST(IngredientsTest, GetterAndSetter){

		typedef LOKI_TYPELIST_2(FeatureBox, FeatureBondset<>) Features;
		typedef ConfigureSystem<VectorInt3,Features,3> Config;
		typedef Ingredients < Config> MyIngredients;
		MyIngredients ingredients;

		EXPECT_EQ(3, ingredients.getMolecules().getMaxConnectivity());

		ingredients.modifyMolecules().setAge(5);

		EXPECT_EQ(5, ingredients.getMolecules().getAge());
		EXPECT_EQ(5, ingredients.modifyMolecules().getAge());

		ingredients.addComment("TestComment\n");
		EXPECT_STREQ("TestComment\n", ingredients.getSumOfComments().c_str());
		ingredients.addComment("TestComment23\n");
		EXPECT_STREQ("TestComment\nTestComment23\n", ingredients.getSumOfComments().c_str());

}

TEST(IngredientsTest, printMetaData){
	typedef LOKI_TYPELIST_2(FeatureBox, FeatureBondset<>) Features;
	typedef ConfigureSystem<VectorInt3,Features,3> Config;
	typedef Ingredients < Config> MyIngredients;

			MyIngredients ingredients;

			stringstream ss;
			ingredients.printMetaData(ss);

			string s;
			ss >> s;

			EXPECT_EQ(ss.str().length(), 92);

}


TEST(IngredientsTest, synchronize_noargument){
	typedef LOKI_TYPELIST_2(FeatureBox, FeatureBondset<>) Features1;
	typedef ConfigureSystem<VectorInt3,Features1,3> Config1;
	typedef Ingredients < Config1> MyIngredients1;

	typedef LOKI_TYPELIST_3(FeatureBox, FeatureBondset<>,FeatureExcludedVolumeSc<>) Features2;
	typedef ConfigureSystem<VectorInt3,Features1,3> Config2;
	typedef Ingredients < Config2> MyIngredients2;


	MyIngredients1 ingredients1;
	MyIngredients2 ingredients2;

	//exception from FeatureBox is exspected, as box size is not set
	EXPECT_THROW(ingredients1.synchronize(),std::runtime_error);
	EXPECT_THROW(ingredients2.synchronize(),std::runtime_error);

	ingredients1.setBoxX(10);
	ingredients1.setBoxY(10);
	ingredients1.setBoxZ(10);
	ingredients1.setPeriodicX(true);
	ingredients1.setPeriodicY(true);

	ingredients2.setBoxX(10);
	ingredients2.setBoxY(10);
	ingredients2.setBoxZ(10);
	ingredients2.setPeriodicX(true);
	ingredients2.setPeriodicY(false);

	//exception from FeatureBox is exspected, as periodicity not fully set
	EXPECT_THROW(ingredients1.synchronize(),std::runtime_error);
	EXPECT_THROW(ingredients2.synchronize(),std::runtime_error);

	ingredients1.setPeriodicZ(false);
	ingredients2.setPeriodicZ(false);

	//now add the bondset and some particles so that effect of synchronize
	//can be checked
	ingredients1.modifyBondset().addBFMclassicBondset();
	ingredients2.modifyBondset().addBFMclassicBondset();

	ingredients1.modifyMolecules().addMonomer(1,1,1);
	ingredients1.modifyMolecules().addMonomer(1,4,6);
	ingredients1.modifyMolecules().addMonomer(7,7,1);
	//this one violates excluded volume, but excluded volume feauture is not used
	ingredients1.modifyMolecules().addMonomer(7,8,2);
	//in ingredients1 this one is ok, but in ingredients2 it will violate box dimensions
	ingredients1.modifyMolecules().addMonomer(3,9,2);

	ingredients2.modifyMolecules().addMonomer(1,1,1);
	ingredients2.modifyMolecules().addMonomer(1,4,6);
	ingredients2.modifyMolecules().addMonomer(7,7,1);
	//this one violates excluded volume
	ingredients2.modifyMolecules().addMonomer(7,8,2);
	//this one violates box dimensions due to periodicity
	ingredients2.modifyMolecules().addMonomer(3,9,2);

	//now make some invalid bonds
	ingredients1.modifyMolecules().connect(0,1);
	ingredients2.modifyMolecules().connect(0,1);


	//again, try to synchronize.
	EXPECT_THROW(ingredients1.synchronize(),std::runtime_error);
	EXPECT_THROW(ingredients2.synchronize(),std::runtime_error);

	//change position of monomer 1, so that bond is valid.
	ingredients1.modifyMolecules()[1].setZ(1);
	ingredients2.modifyMolecules()[1].setZ(1);

	//again, try to synchronize. ingredients1 should work now
	EXPECT_NO_THROW(ingredients1.synchronize());
	EXPECT_THROW(ingredients2.synchronize(),std::runtime_error);

	//now change the conflicting particle coordinate for excluded volume
	//and synchronize again..still fails due to periodicit
	ingredients2.modifyMolecules()[3].setY(5);
	EXPECT_THROW(ingredients2.synchronize(),std::runtime_error);

	//now change the conflicting particle coordinate for periodicity
	//and synchronize again..should work now
	ingredients2.modifyMolecules()[4].setY(5);
	EXPECT_NO_THROW(ingredients2.synchronize());
}

TEST(IngredientsTest, synchronize_withargument){
	typedef LOKI_TYPELIST_2(FeatureBox, FeatureBondset<>) Features1;
	typedef ConfigureSystem<VectorInt3,Features1,3> Config1;
	typedef Ingredients < Config1> MyIngredients1;

	typedef LOKI_TYPELIST_3(FeatureBox, FeatureBondset<>,FeatureExcludedVolumeSc<>) Features2;
	typedef ConfigureSystem<VectorInt3,Features1,3> Config2;
	typedef Ingredients < Config2> MyIngredients2;


	MyIngredients1 ingredients1;
	MyIngredients2 ingredients2;

	//exception from FeatureBox is exspected, as box size is not set
	EXPECT_THROW(ingredients1.synchronize(ingredients1),std::runtime_error);
	EXPECT_THROW(ingredients2.synchronize(ingredients2),std::runtime_error);

	ingredients1.setBoxX(10);
	ingredients1.setBoxY(10);
	ingredients1.setBoxZ(10);
	ingredients1.setPeriodicX(true);
	ingredients1.setPeriodicY(true);

	ingredients2.setBoxX(10);
	ingredients2.setBoxY(10);
	ingredients2.setBoxZ(10);
	ingredients2.setPeriodicX(true);
	ingredients2.setPeriodicY(false);

	//exception from FeatureBox is exspected, as periodicity not fully set
	EXPECT_THROW(ingredients1.synchronize(ingredients1),std::runtime_error);
	EXPECT_THROW(ingredients2.synchronize(ingredients2),std::runtime_error);

	ingredients1.setPeriodicZ(false);
	ingredients2.setPeriodicZ(false);

	//now add the bondset and some particles so that effect of synchronize
	//can be checked
	ingredients1.modifyBondset().addBFMclassicBondset();
	ingredients2.modifyBondset().addBFMclassicBondset();

	ingredients1.modifyMolecules().addMonomer(1,1,1);
	ingredients1.modifyMolecules().addMonomer(1,4,6);
	ingredients1.modifyMolecules().addMonomer(7,7,1);
	//this one violates excluded volume, but excluded volume feauture is not used
	ingredients1.modifyMolecules().addMonomer(7,8,2);
	//in ingredients1 this one is ok, but in ingredients2 it will violate box dimensions
	ingredients1.modifyMolecules().addMonomer(3,9,2);

	ingredients2.modifyMolecules().addMonomer(1,1,1);
	ingredients2.modifyMolecules().addMonomer(1,4,6);
	ingredients2.modifyMolecules().addMonomer(7,7,1);
	//this one violates excluded volume
	ingredients2.modifyMolecules().addMonomer(7,8,2);
	//this one violates box dimensions due to periodicity
	ingredients2.modifyMolecules().addMonomer(3,9,2);

	//now make some invalid bonds
	ingredients1.modifyMolecules().connect(0,1);
	ingredients2.modifyMolecules().connect(0,1);


	//again, try to synchronize.
	EXPECT_THROW(ingredients1.synchronize(ingredients1),std::runtime_error);
	EXPECT_THROW(ingredients2.synchronize(ingredients2),std::runtime_error);

	//change position of monomer 1, so that bond is valid.
	ingredients1.modifyMolecules()[1].setZ(1);
	ingredients2.modifyMolecules()[1].setZ(1);

	//again, try to synchronize. ingredients1 should work now
	EXPECT_NO_THROW(ingredients1.synchronize(ingredients1));
	EXPECT_THROW(ingredients2.synchronize(ingredients2),std::runtime_error);

	//now change the conflicting particle coordinate for excluded volume
	//and synchronize again..still fails due to periodicity
	ingredients2.modifyMolecules()[3].setY(5);
	EXPECT_THROW(ingredients2.synchronize(ingredients2),std::runtime_error);

	//now change the conflicting particle coordinate for periodicity
	//and synchronize again..should work now
	ingredients2.modifyMolecules()[4].setY(5);
	EXPECT_NO_THROW(ingredients2.synchronize(ingredients2));
}
