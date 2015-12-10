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

