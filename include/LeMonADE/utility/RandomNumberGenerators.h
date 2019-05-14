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

#ifndef LEMONADE_UTILITY_RANDOMNUMBERGENERATORS_H
#define LEMONADE_UTILITY_RANDOMNUMBERGENERATORS_H

////////////////////////////////////////////////////////////////////////////////
//disabled by default because of c++11-incompatibility for nvcc in year 2015
//version g++ 4.7.2
//uncomment, if you want to use it
//#define RANDOMNUMBERGENERATOR_ENABLE_CPP11
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include <LeMonADE/utility/R250.h>

/**
 * @file
 * @brief contains class RandomNumberGenerators
 * */

/**
 * @class RandomNumberGenerator
 * @brief Wrapper class for some random number generators,providing consistent interfaces.
 *
 * @details Wrapper class for some random number generators,providing consistent interfaces.
 * The class provides static instances of random number generators, such that
 * the exact same generators can be used everywhere in the project. Note that the
 * generators have to be seeded explicitly, otherwise they will produce always
 * the same numbers. A convenience function seedAll is supplied to randomly
 * seed all supplied generators.
 * Furthermore, the class provides a convenience function for randomly seeding std::rand()
 * from /dev/urandom
 *
 **/

class RandomNumberGenerators
{
	public:
		RandomNumberGenerators();
		~RandomNumberGenerators();

		///////////////////////////////////////////////////////////////
		//functions for getting random numbers

		//R250Engine
		//! returns random unsignet 32 bit integer from R250Engine
		inline uint32_t r250_rand32(){return r250Engine->r250_rand();} //range [0:2e31-1]
		//! returns random double from R250Engine
		inline double r250_drand(){return r250Engine->r250_uniform();} //range [0.0:1.0]

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
		//std::mt19937 (32 bit Mersenne Twister)
		//! //! returns random unsignet 32 bit integer from 32 bit Mersenne Twister
		inline uint32_t mt_rand32(){return mt19937Engine->operator()();} //range[0:2e32-1]
		//! //! returns random double from Mersenne Twister
		inline double mt_drand(){return double(mt_rand32())/double(mt19937Engine->max());}//range [0.0,1.0]
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/


		////////////////////////////////////////////////////////////////
		//functions for seeding the random number generators

		//! methods which take 256+1+1 32 bit seeds to be used as seeds
		void seedAll( uint32_t const * pSeeds );
		void seedAll( std::vector< uint32_t > const & vSeeds );
		void seedAll( std::string const & fname );
		//! randomly initializes all provided RNGs
		void seedAll();
		//R250Engine seeding
		//! randomly seed only R250Engine from /dev/urandom
		void seedR250();
		//! seed R250Engine with array given as argument
		void seedR250( uint32_t const * seedArray );
		//! randomly seed std:rand()
		void seedSTDRAND();
		void seedSTDRAND( uint32_t seed );

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
		//! randomly seed only Mersenne Twister from /dev/urandom
		void seedMT();
		//! seed Mersenne Twister
		void seedMT(uint32_t seed);
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/

	private:

		//instances of random number engines made static so that the
		//exact same RNGs can be used in different places in the code
		//(same seed, state,...)
		//more generators can be placed here, if wanted.

		//! static instance of R250Engine
		static R250* r250Engine;

#ifdef RANDOMNUMBERGENERATOR_ENABLE_CPP11
		static std::mt19937* mt19937Engine;
#endif /*RANDOMNUMBERGENERATOR_ENABLE_CPP11*/
};

#endif /* LEMONADE_UTILITY_RANDOMNUMBERGENERATORS_H */
