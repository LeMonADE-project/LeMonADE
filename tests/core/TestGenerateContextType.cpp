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

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/core/GenerateContextType.h>

#include <iostream>
#include <typeinfo>

using namespace std;

/*****************************************************************************
 * dummy features set up with some dependencies
*****************************************************************************/
    class Feature0;
    class Feature1;
    class Feature2;
    class Feature3;
    class Feature4;
    class Feature5;



    class Feature0:public Feature{
    public:
      typedef LOKI_TYPELIST_1(Feature5) required_features_back;
      template < class IngredientsType > void synchronize(IngredientsType& val){cout<<"0";}
    };

    class Feature1:public Feature{
    public:
      template < class IngredientsType > void synchronize(IngredientsType& val){cout<<"1";}
    };



    class Feature2:public Feature{
    public:
      typedef LOKI_TYPELIST_2(Feature5,Feature1) required_features_front;
      template < class IngredientsType > void synchronize(IngredientsType& val){cout<<"2";}
    };



    class Feature3:public Feature{
    public:
      typedef LOKI_TYPELIST_1(Feature2) required_features_front;
      typedef LOKI_TYPELIST_1(Feature5) required_features_back;
      template < class IngredientsType > void synchronize(IngredientsType& val){cout<<"3";}
    };

    class Feature4:public Feature{
    public:
      typedef LOKI_TYPELIST_1(Feature3) required_features_front;
      template < class IngredientsType > void synchronize(IngredientsType& val){cout<<"4";}
    };

    class Feature5:public Feature{
    public:
      template < class IngredientsType > void synchronize(IngredientsType& val){cout<<"5";}
    };

/*
 * test FullExpand
 * */
TEST(GenerateContexTypeTest, FullExpand)
{
  typedef LOKI_TYPELIST_1(Feature4) Features;
  typedef FullExpand<Features>::Result ExpandedFeatures;
  unsigned int nFeatures= ::Loki::TL::Length<ExpandedFeatures>::value;

  EXPECT_EQ(6,nFeatures);

  typedef LOKI_TYPELIST_6(Feature5,Feature1,Feature2,Feature3,Feature5,Feature4) ExpectedExpanded;
  EXPECT_EQ(typeid(ExpandedFeatures), typeid(ExpectedExpanded));

  EXPECT_EQ(typeid( ::Loki::NullType),typeid(FullExpand< ::Loki::NullType>::Result));
}

/*
 * test InsertBackRequests
 * */
TEST(GenerateContexTypeTest,InsertBackRequests){
  typedef LOKI_TYPELIST_2(Feature0,Feature3) Features;
  typedef InsertBackRequests<Features>::Result ExpandedFeatures;
  typedef LOKI_TYPELIST_4(Feature0,Feature5,Feature3,Feature5) ExpectedExpanded;

  EXPECT_EQ(typeid(ExpandedFeatures),typeid(ExpectedExpanded));
  EXPECT_EQ(typeid( ::Loki::NullType),typeid(InsertBackRequests< ::Loki::NullType>::Result));
}


/*
 * test InsertFeatureRequests
 * */
TEST(GenerateContexTypeTest,InsertFeatureRequests){
  typedef LOKI_TYPELIST_1(Feature4) FeaturesA;
  typedef LOKI_TYPELIST_2(Feature2,Feature5) FeaturesB;
  typedef LOKI_TYPELIST_2(Feature5,Feature0) FeaturesC;
  typedef InsertFeatureRequests<FeaturesA>::Result ExpandedFeaturesA;
  typedef InsertFeatureRequests<FeaturesB>::Result ExpandedFeaturesB;
  typedef InsertFeatureRequests<FeaturesC>::Result ExpandedFeaturesC;

  typedef LOKI_TYPELIST_5(Feature1,Feature2,Feature3,Feature5,Feature4) ExpectedExpandedA;
  typedef LOKI_TYPELIST_3(Feature5,Feature1,Feature2) ExpectedExpandedB;
  typedef LOKI_TYPELIST_2(Feature0,Feature5) ExpectedExpandedC;

  EXPECT_EQ(typeid(ExpandedFeaturesA),typeid(ExpectedExpandedA));
  EXPECT_EQ(typeid(ExpandedFeaturesB),typeid(ExpectedExpandedB));
  EXPECT_EQ(typeid(ExpandedFeaturesC),typeid(ExpectedExpandedC));
}

