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
#include <LeMonADE/core/ConsistencyCheck.h>

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

TEST(ConsistencyCheckTest,Errorcode){
  typedef LOKI_TYPELIST_5(Feature1,Feature2,Feature3,Feature4,Feature5) FeaturesA;
  typedef LOKI_TYPELIST_5(Feature2,Feature3,Feature4,Feature5,Feature1) FeaturesB;
  typedef LOKI_TYPELIST_5(Feature2,Feature4,Feature5,Feature1,Feature3) FeaturesC;

  EXPECT_EQ(1,(ConsistencyCheck<FeaturesA, ::Loki::NullType>::MY_ERRORSTATE) );
  EXPECT_EQ(2,(ConsistencyCheck<FeaturesB, ::Loki::NullType>::MY_ERRORSTATE) );
  EXPECT_EQ(3,(ConsistencyCheck<FeaturesC, ::Loki::NullType>::MY_ERRORSTATE) );
  EXPECT_EQ(0,(ConsistencyCheck< ::Loki::NullType, ::Loki::NullType >::MY_ERRORSTATE) );
}
