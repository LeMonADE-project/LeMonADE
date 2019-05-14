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

#ifndef LEMONADE_UTILITY_R250_H
#define LEMONADE_UTILITY_R250_H

#define R250_RANDOM_PREFETCH			256
#define R250_RAND_NORMALIZE			2.3283064370807974e-10

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <iomanip>

/**
 * @file
 * @brief R250 random number generator
 * */
/**
 * @class R250
 *
 * @brief Random number generator (RNG) for the LeMonADe-project based on Marco's implementation
 *
 * @details Random number generator (RNG) for the LeMonADe-project based on
 * Marco's implementation. The quality of the numbers was tested using the NIST
 * test suite (http://csrc.nist.gov/groups/ST/toolkit/rng/documentation_software.html)
 * and showed good results over a sample of more than 1e9 bits.
 * Reference:Kirkpatrick and Stoll. Journal of Computational Physics, Vol 40, No. 2, 1981.
 * The user has to explicitly seed this random number generator, because it
 * is initially set up in a default state that reproduces a deterministic sequence.
 * For seeding, the function loadRandomState() is supplied.
 */

class R250
{
private:

	//! array holding the 256 random numbers
	uint32_t array[R250_RANDOM_PREFETCH];
	//! points at the next random number from the array
	uint32_t *dice;
	//! always points to end of array.
	uint32_t *arrayEnd;
	//! used in random number generation algorithm
	uint32_t *pos;
	//! used in random number generation algorithm
	uint32_t *other147;
	//! used in random number generation algorithm
	uint32_t *other250;

public:
	//! Constructor
	R250();

	//! returns a random 32bit unsigned integer
	inline uint32_t	r250_rand();
	//! returns a random double in range[0,1]
	inline double r250_uniform();
	//! prints the current numbers in the random number array to std::cout
	void printState();
	//! randomly seed the internal state array from /dev/urandom
	void loadRandomState();
	//! initializes the generator with a predefined state, such that numbers can be reproduced.
	void loadDefaultState();
	//! initializes the internal state array with 256 values from argument stateArray.
	void setState( uint32_t const * stateArray );

private:
	//! applies the random number algorithm to the internal state array (next 256 numbers are generated)
	void refresh();

};


///// definition of inline members /////////////////////////////////
inline uint32_t R250::r250_rand()
{
	if(dice - array == R250_RANDOM_PREFETCH) refresh();
	return *(dice++);
}

inline double R250::r250_uniform()
{
	return ((double)r250_rand())*R250_RAND_NORMALIZE;
}

#endif
