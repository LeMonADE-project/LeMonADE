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

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
#   include <random>
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/
#include <sstream>
#include <stdexcept>
#include <vector>

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

/**
 * Uses seeds from a given array. Currently needs 256+1 or 256+1+1 seeds.
 * Very unsafe. Better to use the std::vector function which checks for enough
 * input!
 */
void RandomNumberGenerators::seedAll( uint32_t const * pSeeds )
{
    uint32_t const * curPos = pSeeds;

	//seed the standard rand
	seedSTDRAND( *( curPos++ ) );
	//randomly initialize all provided RNGs
	seedR250( curPos );
    curPos += R250_RANDOM_PREFETCH;

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
	seedMT( *( curPos++ ) );
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/
}

/**
 * Uses seeds from a given array. Currently needs 256+1 or 256+1+1 seeds.
 * Very unsafe. Better to use the std::vector function which checks for enough
 * input!
 */
void RandomNumberGenerators::seedAll( std::vector< uint32_t > const & vSeeds )
{
    uint32_t const * curPos = &vSeeds[0];
    size_t nSeedsRemaining = vSeeds.size();

	//seed the standard rand
    if ( ! ( nSeedsRemaining >= 1 ) )
        throw std::invalid_argument( "[seedAll] Not enough seeds given" );
    --nSeedsRemaining;
	seedSTDRAND( *( curPos++ ) );


	//randomly initialize all provided RNGs
    if ( ! ( nSeedsRemaining >= R250_RANDOM_PREFETCH ) )
        throw std::invalid_argument( "[seedAll] Not enough seeds given" );
	seedR250( curPos );
    nSeedsRemaining -= R250_RANDOM_PREFETCH;
    curPos          += R250_RANDOM_PREFETCH;

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
    if ( ! ( nSeedsRemaining >= 1 ) )
        throw std::invalid_argument( "[seedAll] Not enough seeds given" );
    --nSeedsRemaining;
	seedMT( *( curPos++ ) );
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/
}

/**
 * Reads random seeds from an input file. The file is supposed to contain raw
 * random data, e.g. created with:
 *  head -c $(( (256+2)*4 )) /dev/urandom > seeds.dat
 */
void RandomNumberGenerators::seedAll( std::string const & fname )
{
    std::ifstream fileSeeds( fname.c_str() );
    if ( fileSeeds.fail() )
    {
        std::stringstream msg;
        msg << "[seedAll] Couldn't open seed file '" << fname << "'";
        throw std::invalid_argument( msg.str() );
    }
    /* seek to end to get file size in bytes */
    fileSeeds.seekg( 0, std::ios_base::end );
    size_t const nBytes = fileSeeds.tellg();
    fileSeeds.seekg( 0 );   // set to beginning again
    if ( fileSeeds.fail() )
    {
        std::stringstream msg;
        msg << "[seedAll] An error occured when trying to get file size for seed file '" << fname << "'";
        throw std::invalid_argument( msg.str() );
    }
    size_t const nLongs = nBytes / sizeof( uint32_t );
    std::cerr << "[seedAll] Seed file contains " << nLongs << " uint32_t data\n";

    std::vector< uint32_t > vSeeds( nLongs );
    fileSeeds.read( (char *) &vSeeds[0], nLongs * sizeof( uint32_t ) );
    seedAll( vSeeds );
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
void RandomNumberGenerators::seedR250( uint32_t const * stateArray )
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

void RandomNumberGenerators::seedSTDRAND( uint32_t const seed )
{
    std::srand( seed );
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
void RandomNumberGenerators::seedMT( uint32_t seed )
{
	mt19937Engine->seed(seed);
}
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/

