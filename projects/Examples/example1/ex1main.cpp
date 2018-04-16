/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Hauke Rabbel)
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

/***********************************************************************
 * This example program demonstrates how to use the provided functions
 * for obtaining pseudo random number series.
 * ********************************************************************/

#include <iostream>
#include <LeMonADE/utility/RandomNumberGenerators.h>

int main(int, char**)
{
	/* ****************************************************************
	* The class RandomNumberGenerators provides an interface for
	* seeding the provided random number generators with a random
	* seed, and obtaining random numbers. The special thing about this
	* class is that all instances use the same random number generator.
	* This means, that if the generator is seeded once in the beginning,
	* it is seeded everywhere in the program with the same seed. This can
	* be convenient if one wants to reproduce the exact same
	* trajectory in a simulation twice, for example for testing
	* purposes. The code below shows how to use this.
	*
	* Generally, one should seed the random number generators
	* somewhere at the beginning of the main function in a fashion
	* similar to what is shown below.
	* ***************************************************************/

	//get an instance of the random number generator class
	RandomNumberGenerators randomNumbers;

	//if we draw "random" numbers now, they will be the same every time
	//every time we start the program:
	std::cout<<"**************************************************\n";
	std::cout<<"EXAMPLE 1: These numbers will always be the same: "
	         <<randomNumbers.r250_rand32()<<" "
		 <<randomNumbers.r250_rand32()<<" "
		 <<randomNumbers.r250_rand32()<<" "
		 <<randomNumbers.r250_rand32()<<std::endl;
	std::cout<<"**************************************************\n";
	//now we can randomly seed all random number generators provided
	//note: this also seeds the std::rand() with a random state
	randomNumbers.seedAll();
	/* ***************************************************************
	* now we can get a series of pseudo random numbers using either
	* std::rand()   (not recommended)   or one of the random number
	* generators provided by the class RandomNumberGenerators. At
	* the moment this is only the R250 generator. It can be accessed
	* in the following ways to get either a 32bit integer or a
	* floating point number (double) in the range [0,1]:
	* **************************************************************/

	int32_t randomInteger=randomNumbers.r250_rand32();
	double randomDouble=randomNumbers.r250_drand();
	std::cout<<"**************************************************\n";
	std::cout<<"EXAMPLE 1: Random numbers. These will be different every time"
	         <<std::endl;
	std::cout<<"randomInteger: "<<randomInteger<<" randomDouble: "
			<<randomDouble<<std::endl<<std::endl;
	std::cout<<"**************************************************\n";
	/* ****************************************************************
	 * If we now create a second instance of RandomNumberGenerators,
	 * the generators will already be seeded, because we seeded
	 * RandomNumberGenerators already above.
	 * ***************************************************************/

	 RandomNumberGenerators randomNumbers2;

	 //now, we print 5 pseudo random numbers to the screen
	 std::cout<<"**************************************************\n";
	 std::cout<<"EXAMPLE 1: numbers from randomNumbers2:\n"
	          <<"These will also be different every time\n";
	 for(int i=0;i<5;i++){
		 std::cout<<randomNumbers2.r250_rand32()<<" ";
	 }
	 std::cout<<"\n**************************************************\n";
	 return 0;
}
