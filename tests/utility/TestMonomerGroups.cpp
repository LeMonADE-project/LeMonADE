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

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/utility/DepthIteratorPredicates.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureBox.h>

#include <list>

/******************************************************************************
 * tests demonstrating the splitting of monomers into different subgroups
 * ****************************************************************************/

using std::list;
using std::vector;

TEST(MonomerGroups, DifferentGroupTypes)
{
  typedef Loki::NullType Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients < Config> MyIngredients;

//   typedef Loki::NullType Features;
//   typedef Molecules<VectorInt3> MyMolecules;
//   typedef Ingredients<MyMolecules,Features> MyIngredients;

  MyIngredients ingredients;

  //create a few connected objects
  ingredients.modifyMolecules().addMonomer(VectorInt3(1,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(2,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(3,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(4,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(5,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(6,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(7,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(8,1,1));
  ingredients.modifyMolecules().addMonomer(VectorInt3(9,1,1));

  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(2,3);
  ingredients.modifyMolecules().connect(3,4);
  ingredients.modifyMolecules().connect(5,6);
  ingredients.modifyMolecules().connect(6,7);
  ingredients.modifyMolecules().connect(7,8);

  //find the groups of connected objects and demonstrate the use of different
  //group containers

  // list - based groups:
  typedef list < int   > Group;
  typedef list < Group > GroupList;

  GroupList connectedObjectsList;
  GroupList linearStrandsList;
  // fills list of groups with connected objects, which belong to a linear strand (including ends)
  fill_connected_groups( ingredients.getMolecules(), linearStrandsList, Group(), belongsToLinearStrand() );

  // vector-based groups:
  typedef vector < MonomerGroup<MyIngredients::molecules_type> > MonomerGroupVector;
  MonomerGroupVector linearStrandsVector;

  fill_connected_groups( ingredients.getMolecules(),
			linearStrandsVector,
			MonomerGroup<MyIngredients::molecules_type>(ingredients.getMolecules()),
			belongsToLinearStrand() );

  //now check that the two groups are the same (subgroups have same size)
  GroupList::iterator li(linearStrandsList.begin());
  MonomerGroupVector::iterator gi(linearStrandsVector.begin());
  while(li!=linearStrandsList.end()){
    EXPECT_EQ(linearStrandsList.size(),linearStrandsVector.size());
    ++li;
    ++gi;
  }
}

TEST(MonomerGroups, CopyGroup)
{
	//define the system
	typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureBox) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;
	//set box
	myIngredients.setBoxX(30);
	myIngredients.setBoxY(30);
	myIngredients.setBoxZ(30);
	//introduce some monomers
	int nMonos=11;
	myIngredients.modifyMolecules().resize(nMonos);

	//set coordinates and connect monomers into chain. also put information
	//on the bonds
	for(int i=0;i<nMonos;i++)
	{
		myIngredients.modifyMolecules()[i].setAllCoordinates(2*i,i,i);
		if(i>0) myIngredients.modifyMolecules().connect(i,i-1,i);
	}

	//also connect monomers 0 and 2
	myIngredients.modifyMolecules().connect(0,2,100);

	//now open a group
	MonomerGroup<Ing::molecules_type> myGroup((myIngredients.getMolecules()));
	//and put some of the monomers into the group
	myGroup.push_back(0);
	myGroup.push_back(1);
	myGroup.push_back(2);
	myGroup.push_back(6);
	myGroup.push_back(8);
	myGroup.push_back(9);
	myGroup.push_back(10);
	//copy the group into a new molecules object
	//this should then contain the correct monomers as well as bonds between
	//them
	Ing::molecules_type mol2=myGroup.copyGroup();

	//now run the tests
	EXPECT_EQ(7,mol2.size());
	EXPECT_EQ(2,mol2.getNumLinks(0));
	EXPECT_EQ(2,mol2.getNumLinks(1));
	EXPECT_EQ(2,mol2.getNumLinks(2));
	EXPECT_EQ(0,mol2.getNumLinks(3)); //previous particle 6
	EXPECT_EQ(1,mol2.getNumLinks(4)); //previous particle 8
	EXPECT_EQ(2,mol2.getNumLinks(5));
	EXPECT_EQ(1,mol2.getNumLinks(6));

	//see if the bond information were also copied
	EXPECT_EQ(1,mol2.getLinkInfo(0,1));
	EXPECT_EQ(2,mol2.getLinkInfo(1,2));
	EXPECT_EQ(100,mol2.getLinkInfo(0,2));
	EXPECT_EQ(9,mol2.getLinkInfo(4,5));
	EXPECT_EQ(10,mol2.getLinkInfo(6,5));



}


TEST(MonomerGroups, ClearGroup)
{
	//define the system
	typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureBox) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;
	//set box
	myIngredients.setBoxX(30);
	myIngredients.setBoxY(30);
	myIngredients.setBoxZ(30);
	//introduce some monomers
	int nMonos=11;
	myIngredients.modifyMolecules().resize(nMonos);

	//set coordinates and connect monomers into chain. also put information
	//on the bonds
	for(int i=0;i<nMonos;i++)
	{
		myIngredients.modifyMolecules()[i].setAllCoordinates(2*i,i,i);
		if(i>0) myIngredients.modifyMolecules().connect(i,i-1,i);
	}

	//also connect monomers 0 and 2
	myIngredients.modifyMolecules().connect(0,2,100);

	//now open a group
	MonomerGroup<Ing::molecules_type> myGroup((myIngredients.getMolecules()));
	//and put some of the monomers into the group
	myGroup.push_back(0);
	myGroup.push_back(1);
	myGroup.push_back(2);
	myGroup.push_back(6);
	myGroup.push_back(8);
	myGroup.push_back(9);
	myGroup.push_back(10);

	EXPECT_EQ(myGroup.size(),7);

	//clear the group
	myGroup.clear();
	EXPECT_EQ(myGroup.size(),0);

	//put some stuff back
	myGroup.push_back(2);
	myGroup.push_back(6);
	myGroup.push_back(8);
	myGroup.push_back(9);
	myGroup.push_back(10);

	EXPECT_EQ(myGroup.size(),5);
	EXPECT_EQ(myGroup.trueIndex(0),2);
	EXPECT_EQ(myGroup.trueIndex(2),8);
	EXPECT_EQ(myGroup.trueIndex(4),10);
}

TEST(MonomerGroups, Erase)
{
	//define the system
	typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureBox) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;
	//set box
	myIngredients.setBoxX(30);
	myIngredients.setBoxY(30);
	myIngredients.setBoxZ(30);
	//introduce some monomers
	int nMonos=11;
	myIngredients.modifyMolecules().resize(nMonos);

	//set coordinates and connect monomers into chain. also put information
	//on the bonds
	for(int i=0;i<nMonos;i++)
	{
		myIngredients.modifyMolecules()[i].setAllCoordinates(2*i,i,i);
		if(i>0) myIngredients.modifyMolecules().connect(i,i-1,i);
	}

	//also connect monomers 0 and 2
	myIngredients.modifyMolecules().connect(0,2,100);

	//now open a group
	MonomerGroup<Ing::molecules_type> myGroup((myIngredients.getMolecules()));
	//and put some of the monomers into the group
	myGroup.push_back(0);
	myGroup.push_back(1);
	myGroup.push_back(2);
	myGroup.push_back(6);
	myGroup.push_back(8);
	myGroup.push_back(9);
	myGroup.push_back(10);

	EXPECT_EQ(myGroup.size(),7);

	//delete stuff from the group
	EXPECT_NO_THROW(myGroup.erase(6));
	EXPECT_THROW(myGroup.erase(7),std::runtime_error);

	EXPECT_EQ(myGroup.size(),6);

	myGroup.erase(2);
	myGroup.erase(2);

	EXPECT_EQ(myGroup.size(),4);

	EXPECT_EQ(myGroup.trueIndex(0),0);
	EXPECT_EQ(myGroup.trueIndex(1),1);
	EXPECT_EQ(myGroup.trueIndex(2),8);
	EXPECT_EQ(myGroup.trueIndex(3),9);



}


TEST(MonomerGroups, RemoveFromGroup)
{
	//define the system
	typedef LOKI_TYPELIST_2(FeatureBondset<>,FeatureBox) Features;
	typedef ConfigureSystem<VectorInt3,Features> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;
	//set box
	myIngredients.setBoxX(30);
	myIngredients.setBoxY(30);
	myIngredients.setBoxZ(30);
	//introduce some monomers
	int nMonos=11;
	myIngredients.modifyMolecules().resize(nMonos);

	//set coordinates and connect monomers into chain. also put information
	//on the bonds
	for(int i=0;i<nMonos;i++)
	{
		myIngredients.modifyMolecules()[i].setAllCoordinates(2*i,i,i);
		if(i>0) myIngredients.modifyMolecules().connect(i,i-1,i);
	}

	//also connect monomers 0 and 2
	myIngredients.modifyMolecules().connect(0,2,100);

	//now open a group
	MonomerGroup<Ing::molecules_type> myGroup((myIngredients.getMolecules()));
	//and put some of the monomers into the group
	myGroup.push_back(0);
	myGroup.push_back(1);
	myGroup.push_back(2);
	myGroup.push_back(6);
	myGroup.push_back(8);
	myGroup.push_back(9);
	myGroup.push_back(10);

	EXPECT_EQ(myGroup.size(),7);

	//delete stuff from the group
	EXPECT_NO_THROW(myGroup.removeFromGroup(6));
	EXPECT_NO_THROW(myGroup.removeFromGroup(10));
	EXPECT_THROW(myGroup.removeFromGroup(7),std::runtime_error);
	EXPECT_THROW(myGroup.removeFromGroup(3),std::runtime_error);

	EXPECT_EQ(myGroup.size(),5);

	EXPECT_EQ(myGroup.trueIndex(0),0);
	EXPECT_EQ(myGroup.trueIndex(1),1);
	EXPECT_EQ(myGroup.trueIndex(2),2);
	EXPECT_EQ(myGroup.trueIndex(3),8);
	EXPECT_EQ(myGroup.trueIndex(4),9);



}