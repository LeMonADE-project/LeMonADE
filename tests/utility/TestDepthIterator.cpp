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

#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/utility/MonomerGroup.h>

#include <LeMonADE/feature/FeatureAttributes.h>

using namespace std;

class TestDepthIterator: public ::testing::Test{
public:

  //redirect cout output
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };

private:
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
};


TEST_F( TestDepthIterator, Iteration )
{
	typedef Molecules < VectorInt3, 3 > GraphType;

	GraphType graph;

	// Three graphs containing elements (1,2,4,5,9), (0,3,6,7), and (8) :
	//
	//  5 - 1   0 - 3 - 6   8
	//  |   |       |
	//  2 - 4       7
	//      |
	//      9

	graph.resize(10);

	graph.connect(5,1);
	graph.connect(1,4);
	graph.connect(4,2);
	graph.connect(2,5);
	graph.connect(4,9);

	graph.connect(0,3);
	graph.connect(3,7);
	graph.connect(3,6);

	set < int > visited;

	GraphIteratorDepthFirst < GraphType > iter1(graph, &visited);

	EXPECT_EQ(iter1.getVertexIdx(),0);

	++iter1; EXPECT_EQ(iter1.getVertexIdx(),3);
	++iter1; EXPECT_EQ(iter1.getVertexIdx(),7);
	++iter1; EXPECT_EQ(iter1.getVertexIdx(),6);
	++iter1; EXPECT_TRUE(iter1.isEnd());

	GraphIteratorDepthFirst < GraphType > iter2(graph, &visited);

	EXPECT_EQ(iter2.getVertexIdx(),1);

	++iter2; EXPECT_EQ(iter2.getVertexIdx(),5);
	++iter2; EXPECT_EQ(iter2.getVertexIdx(),2);
	++iter2; EXPECT_EQ(iter2.getVertexIdx(),4);
	++iter2; EXPECT_EQ(iter2.getVertexIdx(),9);
	++iter2; EXPECT_TRUE(iter2.isEnd());

	GraphIteratorDepthFirst < GraphType > iter3(graph, &visited);

	EXPECT_EQ(iter3.getVertexIdx(),8);

	++iter3; EXPECT_TRUE(iter3.isEnd());

	/// @todo test predicated walks also.
}

TEST_F( TestDepthIterator, DepthIteratorPredicates )
{
  typedef LOKI_TYPELIST_1(FeatureAttributes<>) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> IngredientsType;
  
  IngredientsType ingredients;
  
  // add some tags
  ingredients.modifyMolecules().resize(8);
  ingredients.modifyMolecules()[0].setAttributeTag(1);
  ingredients.modifyMolecules()[1].setAttributeTag(2);
  ingredients.modifyMolecules()[2].setAttributeTag(5);
  ingredients.modifyMolecules()[3].setAttributeTag(4);
  ingredients.modifyMolecules()[4].setAttributeTag(3);
  ingredients.modifyMolecules()[5].setAttributeTag(1);
  ingredients.modifyMolecules()[6].setAttributeTag(3);
  ingredients.modifyMolecules()[7].setAttributeTag(3);
  
  //add some positions
  ingredients.modifyMolecules()[0].setAllCoordinates(0,0,0);
  ingredients.modifyMolecules()[1].setAllCoordinates(0,0,2);
  ingredients.modifyMolecules()[2].setAllCoordinates(0,0,4);
  ingredients.modifyMolecules()[3].setAllCoordinates(0,0,6);
  ingredients.modifyMolecules()[4].setAllCoordinates(0,0,8);
  ingredients.modifyMolecules()[5].setAllCoordinates(0,0,10);
  ingredients.modifyMolecules()[6].setAllCoordinates(0,0,12);
  ingredients.modifyMolecules()[7].setAllCoordinates(0,0,14);
  
  //connect all to stay within one molecule
  ingredients.modifyMolecules().connect(0,1);
  ingredients.modifyMolecules().connect(1,2);
  ingredients.modifyMolecules().connect(2,3);
  ingredients.modifyMolecules().connect(3,4);
  ingredients.modifyMolecules().connect(4,5);
  ingredients.modifyMolecules().connect(5,6);
  ingredients.modifyMolecules().connect(6,7);
  
  // setup monomer groups vector
  typedef std::vector < MonomerGroup<IngredientsType::molecules_type> > MonomerGroupVectorType;
  
  MonomerGroupVectorType groupsVector;
  
  // use fill_connected_groups to create monomer groups 
  // test hasThisType of unconnected monomers
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasThisType<1>());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(1).size(),1);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),5);
  
  // test hasThisType of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasThisType<3>());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),4);
  EXPECT_EQ(groupsVector.at(1).size(),2);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),6);
  EXPECT_EQ(groupsVector.at(1).trueIndex(1),7);
  
  // test hasOneOfTheseTwoTypes of unconnected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasOneOfTheseTwoTypes<2,4>());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),1);
  EXPECT_EQ(groupsVector.at(1).size(),1);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),3);
  
  // test hasOneOfTheseTwoTypes of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasOneOfTheseTwoTypes<1,3>());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(1).size(),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(1),5);
  EXPECT_EQ(groupsVector.at(1).trueIndex(2),6);
  EXPECT_EQ(groupsVector.at(1).trueIndex(3),7);
  
  // test notOfBothTypes of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), notOfType<3>());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(1).size(),1);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),5);
  
  // test notOfBothTypes of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), notOfBothTypes<2,4>());
  EXPECT_EQ(groupsVector.size(),3);
  EXPECT_EQ(groupsVector.at(0).size(),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(1).size(),1);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),2);
  EXPECT_EQ(groupsVector.at(2).size(),4);
  EXPECT_EQ(groupsVector.at(2).trueIndex(0),4);
  EXPECT_EQ(groupsVector.at(2).trueIndex(1),5);
  EXPECT_EQ(groupsVector.at(2).trueIndex(2),6);
  EXPECT_EQ(groupsVector.at(2).trueIndex(3),7);
  
  // test alwaysTrue of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), alwaysTrue());
  EXPECT_EQ(groupsVector.size(),1);
  EXPECT_EQ(groupsVector.at(0).size(),8);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
  EXPECT_EQ(groupsVector.at(0).trueIndex(7),7);
  
  // test belongsToLinearStrand of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), belongsToLinearStrand());
  EXPECT_EQ(groupsVector.size(),1);
  EXPECT_EQ(groupsVector.at(0).size(),8);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
  EXPECT_EQ(groupsVector.at(0).trueIndex(7),7);
  
  // test hasBonds of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasBonds());
  EXPECT_EQ(groupsVector.size(),1);
  EXPECT_EQ(groupsVector.at(0).size(),8);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
  EXPECT_EQ(groupsVector.at(0).trueIndex(7),7);
  
  // change connectivities and redo the last three tests
  ingredients.modifyMolecules().disconnect(6,7);
  
  // test alwaysTrue of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), alwaysTrue());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),7);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
  EXPECT_EQ(groupsVector.at(1).size(),1);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),7);
  
  // test belongsToLinearStrand of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), belongsToLinearStrand());
  EXPECT_EQ(groupsVector.size(),1);
  EXPECT_EQ(groupsVector.at(0).size(),7);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
  
  // test hasBonds of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasBonds());
  EXPECT_EQ(groupsVector.size(),1);
  EXPECT_EQ(groupsVector.at(0).size(),7);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
 
  // change connectivities and redo the last three tests
  ingredients.modifyMolecules().connect(6,7);
  ingredients.modifyMolecules().disconnect(3,4);
  
  // test alwaysTrue of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), alwaysTrue());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(1).size(),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(1),5);
  EXPECT_EQ(groupsVector.at(1).trueIndex(2),6);
  EXPECT_EQ(groupsVector.at(1).trueIndex(3),7);
  
  // test belongsToLinearStrand of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), belongsToLinearStrand());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(1).size(),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(1),5);
  EXPECT_EQ(groupsVector.at(1).trueIndex(2),6);
  EXPECT_EQ(groupsVector.at(1).trueIndex(3),7);
  
  // test hasBonds of connected monomers
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), hasBonds());
  EXPECT_EQ(groupsVector.size(),2);
  EXPECT_EQ(groupsVector.at(0).size(),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(1).size(),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),4);
  EXPECT_EQ(groupsVector.at(1).trueIndex(1),5);
  EXPECT_EQ(groupsVector.at(1).trueIndex(2),6);
  EXPECT_EQ(groupsVector.at(1).trueIndex(3),7);
 
  // add a very simple branch
  ingredients.modifyMolecules().resize(9);
  ingredients.modifyMolecules()[8].setAttributeTag(3);
  ingredients.modifyMolecules()[0].setAllCoordinates(0,2,12);
  ingredients.modifyMolecules().connect(6,8);
  ingredients.modifyMolecules().connect(3,4);
  
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), belongsToLinearStrand());
  EXPECT_EQ(groupsVector.size(),3);
  EXPECT_EQ(groupsVector.at(0).size(),6);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(1).size(),1);
  EXPECT_EQ(groupsVector.at(1).trueIndex(0),7);
  EXPECT_EQ(groupsVector.at(2).size(),1);
  EXPECT_EQ(groupsVector.at(2).trueIndex(0),8);
  
  // now always true and belongsToLinearStrand should differ
  groupsVector.clear();
  fill_connected_groups(ingredients.getMolecules(), groupsVector, ingredients.getMolecules(), alwaysTrue());
  EXPECT_EQ(groupsVector.size(),1);
  EXPECT_EQ(groupsVector.at(0).size(),9);
  EXPECT_EQ(groupsVector.at(0).trueIndex(0),0);
  EXPECT_EQ(groupsVector.at(0).trueIndex(1),1);
  EXPECT_EQ(groupsVector.at(0).trueIndex(2),2);
  EXPECT_EQ(groupsVector.at(0).trueIndex(3),3);
  EXPECT_EQ(groupsVector.at(0).trueIndex(4),4);
  EXPECT_EQ(groupsVector.at(0).trueIndex(5),5);
  EXPECT_EQ(groupsVector.at(0).trueIndex(6),6);
  EXPECT_EQ(groupsVector.at(0).trueIndex(7),7);
  EXPECT_EQ(groupsVector.at(0).trueIndex(8),8);
}

/// @todo test generate_connected_groups.

