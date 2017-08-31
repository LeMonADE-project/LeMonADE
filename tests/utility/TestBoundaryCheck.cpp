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

#include <LeMonADE/utility/BoundaryCheck.h>


TEST(IndexOutOfBoundsException, Throwing )
{
	char array[] = {'a','b','b','d'};
	EXPECT_THROW( throw IndexOutOfBoundsException(array,10), std::out_of_range  );
}

TEST(BoundaryCheckPolicies, Dynamic )
{
  	uint length = 4;
	std::vector<int> array(length,0);

	if ( length > 2 )
	{
		EXPECT_NO_THROW ( CheckDynamicBounds :: isValidIdx (array,length-2) );
		EXPECT_NO_THROW ( CheckDynamicBounds :: isValidIdx (array,length-1) );
	}
	EXPECT_THROW ( CheckDynamicBounds :: isValidIdx (array,length   ), IndexOutOfBoundsException );
	EXPECT_THROW ( CheckDynamicBounds :: isValidIdx (array,length+1 ), IndexOutOfBoundsException );

	array.resize(2*length,0);

	if ( length > 0 )
	{
		EXPECT_NO_THROW ( CheckDynamicBounds :: isValidIdx (array,0) );
	}
	EXPECT_THROW ( CheckDynamicBounds :: isValidIdx (array,2*length ), IndexOutOfBoundsException );
	EXPECT_THROW ( CheckDynamicBounds :: isValidIdx (array,-1 ), IndexOutOfBoundsException );
}

TEST(BoundaryCheckPolicies, Static )
{
  	const uint length = 4;
	int array[length] = {0,0,0,0};

	typedef CheckStaticBounds <length,false> BundaryCheckPolicy ;
	typedef CheckStaticBounds <length,true>  VerboseBundaryCheckPolicy ;

	if ( length > 2 )
	{
		EXPECT_NO_THROW ( BundaryCheckPolicy :: isValidIdx (array,length-2) );
		EXPECT_NO_THROW ( BundaryCheckPolicy :: isValidIdx (array,length-1) );
	}

	EXPECT_THROW ( BundaryCheckPolicy :: isValidIdx (array,length   ), std::out_of_range );
	EXPECT_THROW ( BundaryCheckPolicy :: isValidIdx (array,length+1 ), std::out_of_range );
	EXPECT_THROW ( BundaryCheckPolicy :: isValidIdx (array,-1       ), std::out_of_range );
	EXPECT_THROW ( VerboseBundaryCheckPolicy:: isValidIdx (array, -1 ), std::out_of_range );

}

TEST(BoundaryCheckPolicies, NoCheck )
{
	uint length = 4;
	std::vector<int> array(length,0);
	EXPECT_NO_THROW ( DontCheckBounds :: isValidIdx (array,0) );
	EXPECT_NO_THROW ( DontCheckBounds :: isValidIdx (array,1) );
	EXPECT_NO_THROW ( DontCheckBounds :: isValidIdx (array,length-1) );
	EXPECT_NO_THROW ( DontCheckBounds :: isValidIdx (array,length  ) );
	EXPECT_NO_THROW ( DontCheckBounds :: isValidIdx (array,length+1) );
	EXPECT_NO_THROW ( DontCheckBounds :: isValidIdx (array,-1) );

}

