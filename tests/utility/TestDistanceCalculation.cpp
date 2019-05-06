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
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/utility/DistanceCalculation.h>

using namespace Lemonade;

class TestDistanceCalculation: public ::testing::Test{
public:

  
  //used for testing another tag type 
  typedef LOKI_TYPELIST_1(FeatureBox) Features;
  typedef ConfigureSystem<VectorInt3,Features> Config;
  typedef Ingredients<Config> Ing;
  

  //redirect cout output 
  /*
  virtual void SetUp(){
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  };

  //restore original output
  virtual void TearDown(){
    std::cout.rdbuf(originalBuffer);
  };

  /* un/-comment cout */

  Ing ing;
private:
  
  std::streambuf* originalBuffer;
  std::ostringstream tempStream;
};


TEST_F( TestDistanceCalculation, PowerOfTwoBoxes )
{
  ing.setBoxX(8);
  ing.setPeriodicX(true);
  ing.setBoxY(8);
  ing.setPeriodicY(true);
  ing.setBoxZ(8);
  ing.setPeriodicZ(true);
  ing.modifyMolecules().resize(5);
  ing.modifyMolecules()[0].modifyVector3D().setAllCoordinates(1,1,1);
  ing.modifyMolecules()[1].modifyVector3D().setAllCoordinates(2,1,1);
  ing.modifyMolecules()[2].modifyVector3D().setAllCoordinates(2,11,11); //(2,3,3)
  ing.modifyMolecules()[3].modifyVector3D().setAllCoordinates(66,513,497); //(2,1,1)
  ing.modifyMolecules()[4].modifyVector3D().setAllCoordinates(1,1,-3);
  
  Ing::molecules_type const& mol = ing.getMolecules();
  int dist;
  VectorInt3 vec;
  dist = MinImageDistanceForPowerOfTwo( mol[0],mol[1], ing );
  vec  = MinImageVectorForPowerOfTwo  ( mol[0],mol[1], ing );
  EXPECT_EQ(1,dist);
  EXPECT_EQ(VectorInt3(1,0,0),vec);
  
  dist = MinImageDistanceForPowerOfTwo( mol[0],mol[2], ing );
  vec  = MinImageVectorForPowerOfTwo  ( mol[0],mol[2], ing );
  EXPECT_EQ(3,dist);
  EXPECT_EQ(VectorInt3(1,2,2),vec);
  vec  = MinImageVectorForPowerOfTwo  ( mol[2],mol[0], ing );
  EXPECT_EQ(VectorInt3(-1,-2,-2),vec);  

  dist = MinImageDistanceForPowerOfTwo( mol[0],mol[3], ing );
  vec  = MinImageVectorForPowerOfTwo  ( mol[0],mol[3], ing );
  EXPECT_EQ(1,dist);
  EXPECT_EQ(VectorInt3(1,0,0),vec);
  
  dist = MinImageDistanceForPowerOfTwo( mol[0],mol[4], ing );
  vec  = MinImageVectorForPowerOfTwo  ( mol[0],mol[4], ing );
  EXPECT_EQ(4,dist);
  EXPECT_EQ(VectorInt3(0,0,-4),vec);

  dist = MinImageDistanceForPowerOfTwo( mol[4],mol[0], ing );
  vec  = MinImageVectorForPowerOfTwo  ( mol[4],mol[0], ing );
  EXPECT_EQ(4,dist);
  EXPECT_EQ(VectorInt3(0,0,-4),vec);

  // check exceptions
  ing.setBoxX(9);
  EXPECT_ANY_THROW(MinImageDistanceForPowerOfTwo( mol[0],mol[1], ing ));
  ing.setBoxX(8);
  EXPECT_NO_THROW(MinImageDistanceForPowerOfTwo( mol[0],mol[1], ing ));

  ing.setPeriodicX(false);
  dist = MinImageDistanceForPowerOfTwo( mol[0],mol[1], ing );
  EXPECT_EQ(1,dist);
  vec  = MinImageVectorForPowerOfTwo  ( mol[0],mol[1], ing );
  EXPECT_EQ(VectorInt3(1,0,0),vec);

}

TEST_F( TestDistanceCalculation, NonPowerOfTwoBoxes )
{
  ing.setBoxX(13);
  ing.setPeriodicX(true);
  ing.setBoxY(13);
  ing.setPeriodicY(true);
  ing.setBoxZ(13);
  ing.setPeriodicZ(true);
  ing.modifyMolecules().resize(6);
  ing.modifyMolecules()[0].modifyVector3D().setAllCoordinates(1,1,1);
  ing.modifyMolecules()[1].modifyVector3D().setAllCoordinates(2,1,1);
  ing.modifyMolecules()[2].modifyVector3D().setAllCoordinates(2,11,11); //(-1,3,3)
  ing.modifyMolecules()[3].modifyVector3D().setAllCoordinates(1,1,54); //(0,0,-1)
  ing.modifyMolecules()[4].modifyVector3D().setAllCoordinates(1,1,-3);
  ing.modifyMolecules()[5].modifyVector3D().setAllCoordinates(1,1,-52); //(0,0,1)
  
  Ing::molecules_type const& mol = ing.getMolecules();
  double dist;
  VectorInt3 vec;
  dist = MinImageDistance( mol[0],mol[1], ing );
  EXPECT_DOUBLE_EQ(1,dist);
  vec  = MinImageVector  ( mol[0],mol[1], ing );
  EXPECT_EQ(VectorInt3(1,0,0),vec);
  
  dist = MinImageDistance( mol[0],mol[2], ing );
  vec  = MinImageVector  ( mol[0],mol[2], ing );
  EXPECT_NEAR(4.3588989,dist,0.0000001);
  EXPECT_EQ(VectorInt3(1,-3,-3),vec);
  vec  = MinImageVector  ( mol[2],mol[0], ing );
  EXPECT_EQ(VectorInt3(-1,3,3),vec);  

  dist = MinImageDistance( mol[0],mol[3], ing );
  vec  = MinImageVector  ( mol[0],mol[3], ing );
  EXPECT_DOUBLE_EQ(1,dist);
  EXPECT_EQ(VectorInt3(0,0,1),vec);
  
  dist = MinImageDistance( mol[0],mol[4], ing );
  vec  = MinImageVector  ( mol[0],mol[4], ing );
  EXPECT_DOUBLE_EQ(4,dist);
  EXPECT_EQ(VectorInt3(0,0,-4),vec);

  dist = MinImageDistance( mol[0],mol[5], ing );
  vec  = MinImageVector  ( mol[0],mol[5], ing );
  EXPECT_DOUBLE_EQ(1,dist);
  EXPECT_EQ(VectorInt3(0,0,-1),vec);

  // check nonperiodic boundaries
  ing.setPeriodicX(false);
  dist = MinImageDistance( mol[0],mol[1], ing );
  EXPECT_EQ(1,dist);
  vec  = MinImageVector  ( mol[0],mol[1], ing );
  EXPECT_EQ(VectorInt3(1,0,0),vec);

}

TEST_F( TestDistanceCalculation, DoubleValuesArbitraryBoxes )
{
  ing.setBoxX(13);
  ing.setPeriodicX(true);
  ing.setBoxY(13);
  ing.setPeriodicY(true);
  ing.setBoxZ(13);
  ing.setPeriodicZ(false);
  ing.modifyMolecules().resize(4);
  ing.modifyMolecules()[0].modifyVector3D().setAllCoordinates(0,0,0);
  ing.modifyMolecules()[1].modifyVector3D().setAllCoordinates(3,0,0);
  ing.modifyMolecules()[2].modifyVector3D().setAllCoordinates(0,11,11);
  ing.modifyMolecules()[3].modifyVector3D().setAllCoordinates(2,12,11);

  
  Ing::molecules_type const& mol = ing.getMolecules();
  double dist;
  VectorDouble3 CoM1,CoM2,vec;

  // (1.5,0,0)
  CoM1.setX((mol[0].getX()+mol[1].getX())/2.0);
  CoM1.setY((mol[0].getY()+mol[1].getY())/2.0);
  CoM1.setZ((mol[0].getZ()+mol[1].getZ())/2.0);

  // (1,11.5,11)
  CoM2.setX((mol[2].getX()+mol[3].getX())/2.0);
  CoM2.setY((mol[2].getY()+mol[3].getY())/2.0);
  CoM2.setZ((mol[2].getZ()+mol[3].getZ())/2.0);

  //std::cout << CoM1<< " - " << CoM2 << std::endl;

  // diff = (-0.5,-1.5,11)
  vec  = MinImageVector  ( CoM1,CoM2, ing );
  EXPECT_EQ(VectorDouble3(-0.5,-1.5,11),vec);
  vec  = MinImageVector  ( CoM2,CoM1, ing );
  EXPECT_EQ(VectorDouble3(0.5,1.5,-11),vec);
  dist = MinImageDistance( CoM1,CoM2, ing );
  EXPECT_DOUBLE_EQ(sqrt(0.5*0.5+1.5*1.5+11*11),dist);

  ing.setPeriodicZ(true);
  // diff = (0.5,1.5,2)
  vec  = MinImageVector  ( CoM1,CoM2, ing );
  EXPECT_EQ(VectorDouble3(-0.5,-1.5,-2),vec);
  vec  = MinImageVector  ( CoM2,CoM1, ing );
  EXPECT_EQ(VectorDouble3(0.5,1.5,2),vec);
  dist = MinImageDistance( CoM1,CoM2, ing );
  EXPECT_DOUBLE_EQ(sqrt(0.5*0.5+1.5*1.5+4),dist);

  // diff = (3,0,0)-(1.5,0,0) = (1.5,0,0)
  vec  = MinImageVector  ( CoM1, (VectorDouble3)(mol[1]), ing );
  EXPECT_EQ(VectorDouble3(1.5,0,0),vec);
  vec  = MinImageVector  ( (VectorDouble3)(mol[1]), CoM1, ing );
  EXPECT_EQ(VectorDouble3(-1.5,0,0),vec);
  dist = MinImageDistance( CoM1, (VectorDouble3)(mol[1]), ing );
  EXPECT_DOUBLE_EQ(1.5,dist);

  ing.modifyMolecules()[1].modifyVector3D().setAllCoordinates(11,0,0);
  // diff = (11,0,0)-(1.5,0,0) = (-3.5,0,0)
  vec  = MinImageVector  ( CoM1, (VectorDouble3)(mol[1]), ing );
  EXPECT_EQ(VectorDouble3(-3.5,0,0),vec);
  vec  = MinImageVector  ( (VectorDouble3)(mol[1]), CoM1, ing );
  EXPECT_EQ(VectorDouble3(3.5,0,0),vec);
  dist = MinImageDistance( CoM1, (VectorDouble3)(mol[1]), ing );
  EXPECT_DOUBLE_EQ(3.5,dist);
  
  //test for position apart by more than a box size
  VectorDouble3 vec1(3.5,0.5,0.5);
  VectorDouble3 vec2(-16.,-12.5,31.5);
  vec  = MinImageVector  ( vec1, vec2, ing );
  EXPECT_EQ(VectorDouble3(6.5,0,5.0),vec);
  vec  = MinImageVector  ( vec2, vec1, ing );
  EXPECT_EQ(VectorDouble3(-6.5,0,-5.0),vec);
  dist = MinImageDistance( vec1, vec2, ing );
  EXPECT_NEAR(8.2,dist,0.01);
  
}
