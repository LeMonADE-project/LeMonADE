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

//#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/GenerateMonomerType.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/utility/Vector3D.h>


/******************************************************************************
testing the automated monomer decoration
******************************************************************************/
    //define some example monomer extensions
    class MonoFeature1{
    public:
      MonoFeature1():m1(1){}
      int mono1() const {return 1;}
      int m1;
    };

    class MonoFeature3a{
    public:
      MonoFeature3a():m3a1(1),m3a2(2){}
      int mono3a() const {return 31;}
      int m3a1;
      int m3a2;
    };

    class MonoFeature3b{
    public:
      MonoFeature3b():m3b1(1),m3b2(2),m3b3(3){}
      int mono3b() const {return 32;}
      int m3b1;
      int m3b2;
      int m3b3;
    };

    //define some example features containing monomer_extensions
    class Feature1:public Feature{
    public:
      int feature1(){return 1;}
      typedef LOKI_TYPELIST_1(MonoFeature1) monomer_extensions;
    };



    class Feature2:public Feature{
    public:
      int feature2(){return 2;}
    };



    class Feature3:public Feature{
    public:
      int feature3(){return 3;}
      typedef LOKI_TYPELIST_2(MonoFeature3a,MonoFeature3b) monomer_extensions;
    };
TEST(GenerateMonomerTypeTest, AutomatedDecoration)
{




    //now test GenerateMonomerType for some feature lists

    //fist: all above defined dummy-features
    typedef LOKI_TYPELIST_3(Feature1,Feature2,Feature3)	FeaturesAll;
    typedef GenerateMonomerType<VectorInt3,FeaturesAll>::Result myMonomerAll;
    myMonomerAll monomerAllFeatures;

    EXPECT_EQ(1,monomerAllFeatures.mono1());
    EXPECT_EQ(31,monomerAllFeatures.mono3a());
    EXPECT_EQ(32,monomerAllFeatures.mono3b());

    //now only single sets of features
    typedef LOKI_TYPELIST_2(Feature1,Feature2)	Features12;
    typedef GenerateMonomerType<VectorInt3,Features12>::Result myMonomer12;
    myMonomer12 monomerFeatures12;
    EXPECT_EQ(1,monomerFeatures12.mono1());

    typedef LOKI_TYPELIST_2(Feature2,Feature1)	Features21;
    typedef GenerateMonomerType<VectorInt3,Features21>::Result myMonomer21;
    myMonomer12 monomerFeatures21;
    EXPECT_EQ(1,monomerFeatures21.mono1());

    //the order of the features should not make a difference for the type,
    //also, since Feature2 has no monomer extensions, it should make no difference
    //to the type at all. Since the monomer extensions all have some member
    //variables of different size, this can be tested using sizeof

    EXPECT_EQ(sizeof(myMonomer12),sizeof(myMonomer21));

    typedef LOKI_TYPELIST_1(Feature1)	Features1;
    typedef GenerateMonomerType<VectorInt3,Features1>::Result myMonomer1;
    EXPECT_EQ(sizeof(myMonomer12),sizeof(myMonomer1));

    typedef LOKI_TYPELIST_1(Feature2)	Features2;
    typedef GenerateMonomerType<VectorInt3,Features2>::Result myMonomer2;
    EXPECT_EQ(sizeof(VectorInt3),sizeof(myMonomer2));

    EXPECT_EQ(sizeof(VectorInt3)+6*sizeof(int),sizeof(myMonomerAll));

}
