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

#include <vector>
#include <list>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/utility/GroupByProperty.h>
#include <LeMonADE/utility/MonomerGroup.h>

using namespace std;

TEST(GroupByPropertyTest,CreateGroups)
{
  Molecules<VectorInt3> molecules;

  molecules.resize(20);

  //create a tree of monomers
  molecules.connect(0,1);
  molecules.connect(1,2);
  molecules.connect(1,6);
  molecules.connect(2,3);
  molecules.connect(2,4);
  molecules.connect(2,5);
  molecules.connect(6,7);
  molecules.connect(6,8);
  molecules.connect(6,9);

  //create a second tree of the same structure
  molecules.connect(10,11);
  molecules.connect(11,12);
  molecules.connect(11,16);
  molecules.connect(12,13);
  molecules.connect(12,14);
  molecules.connect(12,15);
  molecules.connect(16,17);
  molecules.connect(16,18);
  molecules.connect(16,19);

  //set positions

  molecules[0].setAllCoordinates(0,0,0);
  molecules[1].setAllCoordinates(2,0,0);
  molecules[2].setAllCoordinates(3,-3,0);
  molecules[3].setAllCoordinates(4,-6,0);
  molecules[4].setAllCoordinates(6,-3,0);
  molecules[5].setAllCoordinates(4,0,0);
  molecules[6].setAllCoordinates(3,3,0);
  molecules[7].setAllCoordinates(4,6,0);
  molecules[8].setAllCoordinates(6,3,0);
  molecules[9].setAllCoordinates(4,0,0);

  molecules[10].setAllCoordinates(0,0,10);
  molecules[11].setAllCoordinates(2,0,10);
  molecules[12].setAllCoordinates(3,-3,10);
  molecules[13].setAllCoordinates(4,-6,10);
  molecules[14].setAllCoordinates(6,-3,10);
  molecules[15].setAllCoordinates(4,0,10);
  molecules[16].setAllCoordinates(3,3,10);
  molecules[17].setAllCoordinates(4,6,10);
  molecules[18].setAllCoordinates(6,3,10);
  molecules[19].setAllCoordinates(4,0,10);

  //now group monomers using group_by_property
  //first group into groups with equal number of bonds
  //here, the base type for a group is simply a list of int,
  //while the groups are stored in a std::vector

  std::vector< list < int32_t > > connectivityGroups;
  group_by_property(molecules,connectivityGroups,list <int32_t>(),NumberOfNeighbors());

  //now sort the monomers into groups with respect to their z-coordinate
  //this time use a vector of MonomerGroup objects for the grouping
  //to demonstrate that this works just as well.

  std::vector<MonomerGroup < Molecules <VectorInt3> > > zCoordinateGroups;
  group_by_property(molecules,zCoordinateGroups,MonomerGroup<Molecules<VectorInt3> >(molecules),ZCoordinate());

  //now check if the monomers are correctly sorted into groups
  //start with the connectivityGroups

  //there should be monomers with 1,3 and 4 bonds, thus 3 groups
  EXPECT_EQ(connectivityGroups.size(),3);
  //in the first group there should fourteen monomers with one connection
  EXPECT_EQ(connectivityGroups[0].size(),14);
  size_t monoIndex=*(connectivityGroups[0].begin());
  EXPECT_EQ(1,molecules.getNumLinks(monoIndex));
  //in the second group there should be two monomers with three connections
  EXPECT_EQ(connectivityGroups[1].size(),2);
  monoIndex=*(connectivityGroups[1].begin());
  EXPECT_EQ(3,molecules.getNumLinks(monoIndex));
  //in the third group there should be four monomers with four connections
  EXPECT_EQ(connectivityGroups[2].size(),4);
  monoIndex=*(connectivityGroups[2].begin());
  EXPECT_EQ(4,molecules.getNumLinks(monoIndex));

  //now check the zCoordinateGroups

  //there should be two different groups, each containing one of the trees
  EXPECT_EQ(zCoordinateGroups.size(),2);
  //the x and y coordinates should be the same, the z coordinate shifted by 10
  for(size_t i=0;i<10;++i){
    EXPECT_EQ(zCoordinateGroups[0][i].getX(),zCoordinateGroups[1][i].getX());
    EXPECT_EQ(zCoordinateGroups[0][i].getY(),zCoordinateGroups[1][i].getY());
    EXPECT_EQ(zCoordinateGroups[0][i].getZ(),zCoordinateGroups[1][i].getZ()-10);
  }
}
