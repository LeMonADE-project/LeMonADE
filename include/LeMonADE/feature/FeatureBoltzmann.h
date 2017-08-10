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

#ifndef LEMONADE_FEATURE_FEATUREBOLTZMANN_H
#define LEMONADE_FEATURE_FEATUREBOLTZMANN_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>


/**
 * @file
 * @class FeatureBoltzmann
 *
 * @brief Feature that checks the acceptance of a move according to the Metropolis-criterion.
 *
 * @details This Feature implements the
 * <a href="http://en.wikipedia.org/wiki/Equation_of_State_Calculations_by_Fast_Computing_Machines">Metropolis</a>
 * -criterion by comparing a uniformly random number &zeta; &isin; [0,1) to the probability \a p.
 * The probability \a p is given by the Boltzmann-factor \a p = exp(-&Delta;U/k<sub>B</sub>T) where
 * &Delta;U is the energetic difference in the move, k<sub>B</sub> is the Boltzmann-constant and T the absolute temperature.
 * If random number
 * * &zeta; < \a p &rarr; move is accepted for applying
 * * &zeta; &ge; \a p &rarr; move is rejected
 *
 * This Feature should by used as \n
 * typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;\n
 * for correct functionality.
 *
 * @todo Rename FeatureBoltzmann into FeatureMetropolis???
 **/
class FeatureBoltzmann:public Feature
{
public:
	//! Default constructor (empty)
	FeatureBoltzmann(){}

	//! Default destructor (empty)
	virtual ~FeatureBoltzmann(){}

	//! Checks if move is allowed by the Metropolis-criterion.
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients,const MoveBase& move);

private:
	//! RNG (Random Number Generator) for random numbers. Needs to be seeded in main() or somewhere appropriate.
	RandomNumberGenerators randomNumbers;

};


/**
 * @details The functions checks the overall probability (energetic difference) in the move against
 * the randomly chosen number using the Metropolis-criterion.
 *
 * @param [in] ingredients A reference to the IngredientsType - mainly the system.
 * @param [in] move General Move.
 * @return if move is allowed (true) or rejected (false).
 */
template<class IngredientsType>
bool FeatureBoltzmann::checkMove(const IngredientsType& ingredients, const MoveBase& move)
{
	return( (randomNumbers.r250_drand() < move.getProbability() ) ? true : false);
}


#endif
