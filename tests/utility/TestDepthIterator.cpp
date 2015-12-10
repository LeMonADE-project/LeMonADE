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
#include <LeMonADE/core/Molecules.h>

using namespace std;

TEST( TestGraphIteratorDepthFirst, Iteration )
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

/// @todo test generate_connected_groups.

