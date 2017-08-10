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

#include <LeMonADE/utility/RandomNumberGenerators.h>

R250* RandomNumberGenerators::r250Engine=0;

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
std::mt19937* RandomNumberGenerator::mt19937Engine=0;
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/


/**
 * @brief Standard constructor, sets up the RNG engines and seeds
 * the RNG with fixed (not random numbers)
 *
 * @details The RNG should be seeded with random numbers using seedR250() or seedMT()
 * 			otherwise the RNG is initialized with default values.
 *
 */
RandomNumberGenerators::RandomNumberGenerators()
{

	if(r250Engine==0)
		{
		r250Engine=new R250;
		r250Engine->loadDefaultState();
		}
#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
	if(mt19937Engine==0)
	{
		mt19937Engine=new std::mt19937;
		seedMT(1);
	}
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/

}


RandomNumberGenerators::~RandomNumberGenerators()
{

}

void RandomNumberGenerators::seedAll()
{
	//seed the standard rand
	seedSTDRAND();
	//randomly initialize all provided RNGs
	seedR250();

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
	seedMT();
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/

}

//randomly seed R250Engine
//R250Engine handles getting a random state internally
void RandomNumberGenerators::seedR250()
{
	r250Engine->loadRandomState();
}

//seed R250Engine with array of 256 values
//state array has to contain 256 ideally random values
void RandomNumberGenerators::seedR250(uint32_t* stateArray)
{
	r250Engine->setState(stateArray);
}

//convenience function for randomly seeding std::rand() from /dev/urandom
void RandomNumberGenerators::seedSTDRAND()
{
	std::ifstream urandom("/dev/urandom", std::ios::binary);
	if (urandom.is_open())
	{
		uint32_t* seed=new uint32_t;
		std::ifstream in("/dev/urandom");
		in.read((char*)seed,sizeof(seed));
		std::srand(*seed);
		std::cout<<"RandomNumberGenerators: random seed for std::rand() picked from /dev/urandom:"<<seed<<std::endl;
		urandom.close();
	}
	else
	{
		std::stringstream errormessage;
		errormessage<<"could not generate random seed for std::rand() from urandom..exiting\n";
		throw std::runtime_error(errormessage.str());

	}

}


#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11

//randomly seed Mersenne Twister from /dev/urandom
void RandomNumberGenerators::seedMT()
{
	std::ifstream urandom("/dev/urandom", std::ios::binary);
	if (urandom.is_open())
	{
		uint32_t* seed=new uint32_t;
		std::ifstream in("/dev/urandom");
		in.read((char*)seed,sizeof(seed));
		std::srand(*seed);
		std::cout<<"RandomNumberGenerators: random seed for Mersenne Twister picked from /dev/urandom:"<<seed<<std::endl;
		urandom.close();
	}
	else
	{
		std::stringstream errormessage;
		errormessage<<"could not generate random seed for Mersenne Twister from urandom..exiting\n";
		throw std::runtime_error(errormessage.str());

	}
}

//seed Mersenne Twister
void RandomNumberGenerators::seedMT(uint32_t seed)
{
	mt19937Engine->seed(seed);
}
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/

