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
 * TestConnected.cpp
 *
 *  Created on: 22.04.2013
 *      Author: christoph
 */

#include "gtest/gtest.h"

#include <sstream>
#include <exception>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/core/ConnectedDecorator.h>


TEST(ConnectedTest, GetterAndSetter) {

	Connected<VectorInt3, 5> connected;

	EXPECT_EQ(5, connected.getMaxConnectivity());
	EXPECT_EQ(0, connected.getNumLinks());

	EXPECT_THROW(connected.getNeighborIdx(0), std::runtime_error);

	connected.connect(2);

	EXPECT_EQ(1, connected.getNumLinks());
	EXPECT_EQ(2, connected.getNeighborIdx(0));
}

TEST(ConnectedTest, Operators) {

	Connected<VectorInt3, 5> connected5;
	connected5.connect(2);
	EXPECT_EQ(1, connected5.getNumLinks());

	Connected<VectorInt3, 0> connected0;
	EXPECT_EQ(0, connected0.getNumLinks());

	std::stringstream ss0;
	ss0 << connected0;
	EXPECT_EQ("0 0 0", ss0.str());

	Connected<VectorInt3, 1> connected1;
	connected1.connect(2);
	EXPECT_EQ(1, connected1.getNumLinks());


	Connected<VectorInt3, 2> connected2;
	connected5.connect(100);
	connected5.connect(200);
	connected5.connect(400);
	EXPECT_ANY_THROW(connected2 = connected5);

	EXPECT_THROW((Connected<VectorInt3,0>(connected5)), std::runtime_error);

	EXPECT_THROW((Connected<VectorInt3,1>(connected5)), std::runtime_error);

	EXPECT_THROW((Connected<VectorInt3,3>(connected5)), std::runtime_error);

}

TEST(ConnectedTest, ConnectAndDisconnect) {
	Connected<VectorInt3, 5> connected1;

	connected1.connect(1);
	connected1.connect(2);
	//connected1.connect(2);
	connected1.connect(3);

	EXPECT_EQ(3, connected1.getNumLinks());

	for (int i = 0; i < connected1.getNumLinks(); i++) {
		EXPECT_TRUE(
				1 == connected1.getNeighborIdx(i) ||
				2 == connected1.getNeighborIdx(i) ||
				3 == connected1.getNeighborIdx(i));
	}

	Connected<VectorInt3, 5> connected2(connected1);

	EXPECT_EQ(3, connected2.getNumLinks());

	for (int i = 0; i < connected2.getNumLinks(); i++) {
		EXPECT_TRUE(
				1 == connected1.getNeighborIdx(i) ||
				2 == connected1.getNeighborIdx(i) ||
				3 == connected1.getNeighborIdx(i));
	}

	connected2.disconnect(2);

	EXPECT_EQ(2, connected2.getNumLinks());

	for (int i = 0; i < connected2.getNumLinks(); i++) {
		EXPECT_TRUE( 2 != connected2.getNeighborIdx(i));
	}

	Connected<VectorInt3, 2> connected3;
	EXPECT_ANY_THROW(connected3.connect(-2))<< "negative index is not allowed";
}

//Test for the template specialization for max_connectivity=0
TEST(ConnectedTest, ConnectivityZero){
  //test standard constructor
  Connected<VectorInt3,0> connected0;

  EXPECT_EQ(0,connected0.getMaxConnectivity());
  EXPECT_EQ(0,connected0.getNumLinks());
  EXPECT_THROW(connected0.getNeighborIdx(0),std::runtime_error);

  //est connect and disconnect
  EXPECT_THROW(connected0.connect(1),std::runtime_error);
  //check if everything is still the same
  EXPECT_EQ(0,connected0.getNumLinks());
  EXPECT_THROW(connected0.getNeighborIdx(0),std::runtime_error);
  EXPECT_THROW(connected0.disconnect(1), std::runtime_error);

  //test constructor from vertex
  VectorInt3 pos;
  pos.setAllCoordinates(1,2,3);
  Connected<VectorInt3,0> connected1(pos);

  EXPECT_EQ(1,connected1.getX());
  EXPECT_EQ(2,connected1.getY());
  EXPECT_EQ(3,connected1.getZ());
  EXPECT_EQ(0,connected1.getMaxConnectivity());

  //test conversion constructor
  Connected<VectorInt3,2> connected2;
  connected2.setAllCoordinates(4,5,6);
  Connected<VectorInt3,0> connected3(connected2);
  EXPECT_EQ(4,connected3.getX());
  EXPECT_EQ(5,connected3.getY());
  EXPECT_EQ(6,connected3.getZ());
  connected2.connect(1);
  EXPECT_THROW( (Connected<VectorInt3,0>(connected2)) ,std::runtime_error);

  //test if zero memory is consumed by decorator
  EXPECT_EQ(sizeof(pos),sizeof(connected1));
}

//Test for the template specialization for max_connectivity=1
TEST(ConnectedTest, ConnectivityOne){

  //test standard constructor
  Connected<VectorInt3,1> connected0;
  EXPECT_EQ(1,connected0.getMaxConnectivity());
  EXPECT_EQ(0,connected0.getNumLinks());

  //test connect and disconnect
  connected0.connect(10);
  EXPECT_EQ(1,connected0.getNumLinks());
  EXPECT_EQ(10,connected0.getNeighborIdx(0));

  EXPECT_THROW(connected0.connect(11),std::runtime_error);
  EXPECT_EQ(1,connected0.getNumLinks());

  EXPECT_THROW(connected0.disconnect(11),std::runtime_error);
  EXPECT_EQ(1,connected0.getNumLinks());

  connected0.disconnect(10);
  EXPECT_EQ(0,connected0.getNumLinks());

  //test construction from vertex
  VectorInt3 pos;
  pos.setAllCoordinates(1,2,3);
  Connected<VectorInt3,1> connected1(pos);

  EXPECT_EQ(1,connected1.getX());
  EXPECT_EQ(2,connected1.getY());
  EXPECT_EQ(3,connected1.getZ());
  EXPECT_EQ(1,connected1.getMaxConnectivity());

  //test conversion constructor
  //prepare some objects
  Connected<VectorInt3,3> connected2;
  Connected<VectorInt3,3> connected3;
  connected2.connect(10);
  connected3.connect(20);
  connected3.connect(21);

  connected2.setAllCoordinates(1,2,3);
  //now copy these objects into a new one
  //this one should work
  Connected<VectorInt3,1> connected4(connected2);
  EXPECT_EQ(connected2.getNumLinks(),connected4.getNumLinks());
  EXPECT_EQ(connected2.getNeighborIdx(0),connected4.getNeighborIdx(0));
  EXPECT_EQ(connected2.getX(),connected4.getX());
  EXPECT_EQ(connected2.getY(),connected4.getY());
  EXPECT_EQ(connected2.getZ(),connected4.getZ());
  //this one should not work
  EXPECT_THROW( (Connected<VectorInt3,1>(connected3)),std::runtime_error);


}

