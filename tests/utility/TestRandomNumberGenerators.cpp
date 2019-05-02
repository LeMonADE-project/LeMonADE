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

#include <cstdlib>
#include <cmath>
#include <set>
#include <vector>

#include <LeMonADE/utility/RandomNumberGenerators.h>

//the tests contained here test the interfaces, but also to some very limited
//extend the random numbers. thus, some tests are based on statistics and
//expectation values, and may fail, even if everything works well. if one of
//the tests fails, it may be a good idea to run the test again.


//produces a seed array for r250 and tests for reproducible random numbers
TEST(RandomNumberGeneratorsTest,R250Seeding){

	std::vector<double> numbersDouble;
	std::vector<uint32_t> numbersInt;

	RandomNumberGenerators rng;

	rng.seedSTDRAND();

	uint32_t seedArray[256];
	for(size_t n=0;n<256;n++){
		seedArray[n]=std::rand();
	}

	rng.seedR250(seedArray);

	for(size_t i=0;i<10000;i++){
		numbersDouble.push_back(rng.r250_drand());
		numbersInt.push_back(rng.r250_rand32());
	}

	rng.seedR250(seedArray);

	for(size_t i=0;i<10000;i++){
		EXPECT_EQ(numbersDouble[i],rng.r250_drand());
		EXPECT_EQ(numbersInt[i],rng.r250_rand32());
	}

}

//shows that two instances of RandomNumberGenerators draw from the same
//R250

TEST(RandomNumberGeneratorsTest,R250MultipleInstances){

	std::vector<uint32_t> numbersInt;

	RandomNumberGenerators rng1;
	RandomNumberGenerators rng2;


	rng1.seedSTDRAND();

	uint32_t seedArray[256];
	for(size_t n=0;n<256;n++){
		seedArray[n]=std::rand();
	}

	//seed using rng1
	rng1.seedR250(seedArray);

	//fill the array using rng1 and rng2
	for(size_t i=0;i<5000;i++){
		numbersInt.push_back(rng1.r250_rand32());
		numbersInt.push_back(rng2.r250_rand32());
	}

	//now set back the generator, seed using rng2
	rng2.seedR250(seedArray);

	//compare the random numbers in numbersInt with a new series drawn only
	//from rng1...they should be equal
	for(size_t i=0;i<10000;i++){
		EXPECT_EQ(numbersInt[i],rng1.r250_rand32());
	}
}


//very simple test to see if on average all bits are being used equally
TEST(RandomNumberGeneratorsTest,R250BitsUsage){

	// a statistical test may fail even if everything is fine.
	// to minimize failure rate and keep runtime justifiable the test is repeated if it fails
	// if test fails multiple times something might be seriously wrong and the test does not pass
	uint32_t successCounter(0);
	uint32_t maxFails(10);
	bool testPassed(false);
	
	std::vector<uint32_t> bitUsage(32);

	RandomNumberGenerators rng;
	
	// repeat test is a rate event of failing occurs, but not arbitrarily often
	while( (successCounter<maxFails) && !testPassed){

		//seed using rng
		rng.seedAll();

		//draw 20000 numbers. then, each bit should have been used on average
		//10000 times with a standard deviation of 100. in this loop we count 
		//for each bit the number of times it is used within 20000 random numbers
		for(size_t i=0;i<20000;i++){
			uint32_t randomNumber=rng.r250_rand32();
			for(uint32_t n=1;n<=32;n++){			
				bitUsage[n-1]+=(randomNumber & (uint32_t(std::pow(2,n))-1))>>(n-1);
			}
		}

		std::cout<<"THIS TEST MAY FAIL EVEN IF EVERYTHING IS WORKING PROPERLY\n"
		<<"THE EXPECTED VALUES ARE 10000+-300, WHERE 300 IS 2SIGMA. THAT IS 99.7 PERCENT OF\n"
		<<"THE CASES THE TEST WILL SUCCEED IF EVERYTHING IS OK....IF IT FAILS, TRY RUNNING AGAIN\n";
		//check that all bits have been used are within 3sigma of the expectation
		//value. this test may of course fail

		// store test result in a dummy bool
		bool testBitwisePassed(true);
		for(size_t n=0;n<32;n++){
			//EXPECT_LE(bitUsage[n],(10000+300));
			if( !(bitUsage[n] <= (10000+300)) ) testBitwisePassed=false;
			//EXPECT_GE(bitUsage[n],(10000-300));
			if( !(bitUsage[n] >= (10000-300)) ) testBitwisePassed=false;
		}

		if(testBitwisePassed) testPassed=true;
		successCounter++;
	}
	// use the test routines here
	std::cout << "random number test was executed "<<successCounter<<" times."<<std::endl;
	EXPECT_LT(successCounter,maxFails);
	EXPECT_TRUE(testPassed);
	
}

//this test is somewhat random, it can in principle fail, even if everything is
//working correctly. the probability for the test to fail is roughly 1%.
TEST(RandomNumberGeneratorsTest,R250RangeOfNumbers){
	
	RandomNumberGenerators rng;

	rng.seedAll();

	std::set<double> numbersDouble;
	std::set<uint32_t> numbersInt;

	for(int n=0;n<10000;n++)
	{
		double d=rng.r250_drand();

		EXPECT_TRUE(d>=0.0);
		EXPECT_TRUE(d<=1.0);

		numbersDouble.insert(rng.r250_drand());

		//not checking the range on this one, as it is supposed to fill
		//the whole range of uint32_t anyways. this was checked by NIST
		//test suite
		uint32_t i=rng.r250_rand32();
		numbersInt.insert(i);

	}

	std::cout<<"THIS TEST MAY FAIL EVEN IF EVERYTHING IS WORKING PROPERLY\n"
	<<"THE EXPECTED VALUES ARE 10000, BECAUSE 10000 RANDOM NUMBERS WERE DRAWN\n"
	<<"THE TEST CHECKS IF ALL ARE DIFFERENT, WHICH THEY ARE WITH A PROBABILITY\n"
	<<"OF 98 PERCENT. IF THE TEST FAILS, ONE SHOULD MAYBE RUN IT AGAIN A FEW MORE TIMES\n";
	//check the size of the set. if any two drawn numbers are equal, these
	//tests fail. one would not expect this, but of course it can happen
	//in principle, even if the RNGs work correctly.
	EXPECT_GE(numbersInt.size(),(10000-1));
	EXPECT_GE(numbersDouble.size(),(10000-1));

}