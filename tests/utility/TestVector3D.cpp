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
 * TestVector3D.cpp
 *
 *  Created on: 26.04.2013
 *      Author: christoph
 */

#include "gtest/gtest.h"

#include <sstream>
#include <limits>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/NumericResultTypes.h>

TEST(Vector3DTest, getterAndSetter) {
	VectorInt3 vec;
	EXPECT_EQ(0, vec.getX());
	EXPECT_EQ(0, vec.getY());
	EXPECT_EQ(0, vec.getZ());

	vec.setX(1);
	vec.setY(2);
	vec.setZ(3);

	EXPECT_EQ(1, vec.getX());
	EXPECT_EQ(2, vec.getY());
	EXPECT_EQ(3, vec.getZ());

//	EXPECT_ANY_THROW(vec.SetX(1.1));
//	EXPECT_ANY_THROW(vec.SetY(1.2));
//	EXPECT_ANY_THROW(vec.SetZ(1.3));
	// or compiler warning

	vec.setCoordinate(0, 5);
	EXPECT_EQ(5, vec.getX());
	EXPECT_EQ(5, vec.getCoordinate(0));
	EXPECT_EQ(5, vec[0]);

	vec[0] = 3;
	EXPECT_EQ(3, vec[0]);

	VectorInt3 vec2A(4,9,-3);
	VectorInt3 vec2B(4,9,-3);
	VectorInt3 vecC(-1,2,-4);

	EXPECT_EQ( vec2B.getVector3D(), vec2A.getVector3D());
	EXPECT_EQ( 4, vec2A.getVector3D().getX());
	EXPECT_EQ( 9, vec2A.getVector3D().getY());
	EXPECT_EQ(-3, vec2A.getVector3D().getZ());
	EXPECT_EQ(26, vec2A.getVector3D()*vecC.getVector3D());
	
	EXPECT_NO_THROW(vec2A.modifyVector3D().setX(34));
	EXPECT_EQ( 34, vec2A.getVector3D().getX());

#ifdef DEBUG
	EXPECT_ANY_THROW(vec[5] = int(4));
	EXPECT_ANY_THROW(vec.setCoordinate(5,5));
#endif

}

TEST(Vector3DTest, Constructors) {
	Vector3D<double> vec;
	EXPECT_DOUBLE_EQ(0.0, vec.getX());
	EXPECT_DOUBLE_EQ(0.0, vec.getY());
	EXPECT_DOUBLE_EQ(0.0, vec.getZ());

// 	Vector3D<int> vecI(1);
// 	EXPECT_EQ(1, vecI.getX());
// 	EXPECT_EQ(1, vecI.getY());
// 	EXPECT_EQ(1, vecI.getZ());
//
// 	Vector3D<double> vecD(2.2);
//
// 	Vector3D<int> vecId(vecD);
// 	EXPECT_EQ(2, vecId.getX());
// 	EXPECT_EQ(2, vecId.getY());
// 	EXPECT_EQ(2, vecId.getZ());

// 	Triple<double> pos(1.1, 2.2, 3.3);
// 	Vector3D<int> vecIpos(pos);
// 	EXPECT_EQ(1, vecIpos.getX());
// 	EXPECT_EQ(2, vecIpos.getY());
// 	EXPECT_EQ(3, vecIpos.getZ());

// 	Vector3D<int> vec3Crazy(short(2), double(4.4), bool(true));
// 	EXPECT_EQ(2, vec3Crazy.getX());
// 	EXPECT_EQ(4, vec3Crazy.getY());
// 	EXPECT_EQ(1, vec3Crazy.getZ());

}

TEST(Vector3DTest, Operators) {


	Vector3D<int> vecAssign;
	vecAssign.setAllCoordinates(2, 3, 4);
	EXPECT_EQ(2, vecAssign.getX());
	EXPECT_EQ(3, vecAssign.getY());
	EXPECT_EQ(4, vecAssign.getZ());



//    	Vector3D<int> vec3Crazy(short(2), double(4.4), bool(true));
 	VectorInt3 vec3Crazy(11,22,33);

	VectorFloat3 vecFloat = vec3Crazy;
	EXPECT_DOUBLE_EQ(11.0, vecFloat.getX());
	EXPECT_DOUBLE_EQ(22.0, vecFloat.getY());
	EXPECT_DOUBLE_EQ(33.0, vecFloat.getZ());

 	VectorInt3 vecInt;
	vecInt.setAllCoordinates(2, 3, 4);
	EXPECT_EQ(2, vecInt.getX());
	EXPECT_EQ(3, vecInt.getY());
	EXPECT_EQ(4, vecInt.getZ());

	VectorDouble3 vecDouble(1.1, 2.2, int(4) );

	EXPECT_DOUBLE_EQ(27.3, vecInt * vecDouble + 2.5 );

	EXPECT_DOUBLE_EQ(39.7, vecInt*(1.5*vecDouble) + 2.5 );

	VectorFloat3 offset(0.5f,0.5f,0.5f);

	EXPECT_DOUBLE_EQ(21.5, (vecInt - offset)*(vecInt - offset) + offset*offset);

	EXPECT_FLOAT_EQ(27.3, vecInt * VectorFloat3 ( float ( vecDouble.getX() ), float ( vecDouble.getY() ), float ( vecDouble.getZ() ) ) + 2.5 );

	EXPECT_EQ( typeid ( NumericResultTypes < VectorInt3   , VectorDouble3 > ::product_type ), typeid( double ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorFloat3 , VectorDouble3 > ::product_type ), typeid( double ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorLong3  , VectorFloat3  > ::product_type ), typeid( float ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorFloat3 , VectorLong3   > ::product_type ), typeid( float ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorUint3  , VectorInt3    > ::product_type ), typeid( int ) );
 	EXPECT_EQ( typeid ( NumericResultTypes < VectorUlong3 , VectorInt3    > ::product_type ), typeid( int64_t ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorUlong3 , VectorChar3    > ::product_type ), typeid( int64_t ) );
 	EXPECT_EQ( typeid ( NumericResultTypes < VectorUlong3 , VectorUchar3    > ::product_type ), typeid( uint64_t ) );

	EXPECT_EQ( typeid ( NumericResultTypes < VectorInt3   , VectorDouble3 > ::stronger_type ), typeid( VectorDouble3 ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorFloat3 , VectorDouble3 > ::stronger_type ), typeid( VectorDouble3 ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorLong3  , VectorFloat3  > ::stronger_type ), typeid( VectorFloat3 ) );
	EXPECT_EQ( typeid ( NumericResultTypes < VectorFloat3 , VectorLong3   > ::stronger_type ), typeid( VectorFloat3 ) );
 	EXPECT_EQ( typeid ( NumericResultTypes < VectorUint3  , VectorInt3    > ::stronger_type ), typeid( VectorInt3 ) );
  	EXPECT_EQ( typeid ( NumericResultTypes < VectorUlong3 , VectorInt3    > ::stronger_type ), typeid( VectorLong3 ) );
 	EXPECT_EQ( typeid ( NumericResultTypes < VectorChar3 ,  VectorUlong3  > ::stronger_type ), typeid( VectorLong3 ) );
  	EXPECT_EQ( typeid ( NumericResultTypes < VectorUshort3 ,VectorFloat3  > ::stronger_type ), typeid( VectorFloat3 ) );

}

TEST(Vector3DTest, Normalize) {
	VectorDouble3 a(1.0,2.0,3.0);
	VectorDouble3 b(0.0,0.0,0.0);
	VectorDouble3 c(std::numeric_limits<double>::infinity(),0.0,0.0);
	VectorDouble3 d(-1.0*std::numeric_limits<double>::infinity(),0.0,0.0);

	EXPECT_NO_THROW(a.normalize());
	EXPECT_EQ(a.getLength(),1.0);

	EXPECT_THROW(b.normalize(),std::runtime_error);
	EXPECT_THROW(c.normalize(),std::runtime_error);
	EXPECT_THROW(d.normalize(),std::runtime_error);
}