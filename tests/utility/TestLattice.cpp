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
 * @class Lattice
 * @author Toni Mueller
 * @brief class provides a general lattice 
 */

#include "gtest/gtest.h"
#include <cstdio>
#include <sstream>
#include <LeMonADE/utility/Lattice.h>


class TestLattice: public ::testing::Test{
public:
//   typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes<>, FeatureExcludedVolumeSc<FeatureLatticePowerOfTwo<bool> >) Features;
//   typedef ConfigureSystem<VectorInt3,Features> Config;
//   typedef Ingredients<Config> IngredientsType;
// 
//   IngredientsType ingredients;
//   const IngredientsType& getIngredients() const {return ingredients;}

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
TEST_F(TestLattice, Constructor)
{
  Lattice<> MyLattice;
  MyLattice.setupLattice(4,4,4);
  MyLattice.setLatticeEntry(3,7,12,17);
  Lattice<> MyLattice2(MyLattice);
  MyLattice2.setupLattice();
  EXPECT_EQ(MyLattice2.getLatticeEntry(3,7,12),0);
}
TEST_F(TestLattice, Functionality)
{

  Lattice<> MyLattice;
  MyLattice.setupLattice(4,4,4);
  EXPECT_EQ(MyLattice.getLatticeEntry(0,2,0),0);
  EXPECT_EQ(MyLattice.getLatticeEntry(0,6,0),0);
  MyLattice.setLatticeEntry(0,2,0,17);
  EXPECT_EQ(MyLattice.getLatticeEntry(0,2,0),17);
  EXPECT_EQ(MyLattice.getLatticeEntry(0,6,0),17);
  VectorInt3 Pos(0,2,0);
  EXPECT_EQ(MyLattice.getLatticeEntry(Pos),17);
  VectorInt3 NewPos(12,23,19);
  MyLattice.moveOnLattice(Pos,NewPos);
  EXPECT_EQ(MyLattice.getLatticeEntry(NewPos),17);
  MyLattice.clearLattice();
  EXPECT_EQ(MyLattice.getLatticeEntry(NewPos),0);
 
}
