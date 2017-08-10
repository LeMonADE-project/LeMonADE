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

#ifndef LEMONADE_UPDATER_MOVES_MOVEBASE_H
#define LEMONADE_UPDATER_MOVES_MOVEBASE_H

/**
 * @file
 * @brief contains class Move
 * */
/*****************************************************************************/
/**
 * @class MoveBase
 *
 * @brief Base class for any types of moves inside the simulation.
 *
 * @details All moves must be derived from this base class in order to be processed
 * correctly by the Feature. It is highly recommended (even though not strictly
 * necessary) that all derived moves implement at least the following three
 * member function templates:
 *
 * 	template<class IngredientsType> bool check(IngredientsType& ingredients) const
 * 	{
 * 	return ingredients.checkMove(ingredients,*this);
 * 	}
 *
 *  template<class IngredientsType> void apply(IngredientsType& ingredients)
 * 	{
 * 	ingredients.applyMove(ingredients,*this);
 * 	}
 *
 * 	template <class IngredientsType> void init(const IngredientsType& ingredients){};
 *
 * These three function templates are meant to be called by the user, or by
 * updaters such as UpdaterSimpleSimuator when applying the Monte Carlo move to
 * the simulation system, or checking if they may be applied under the given boundary conditions.
 * Of course, the implementation in specialized moves derived from MoveBase may
 * contain more code than the above minimal examples.
 * For an implementation examples have a look at the classes MoveAddMonomerSc
 * or MoveLocalSc . Examples for theusage of these functions can be found in
 * the classes UpdaterAbstractCreate or UpdaterSimpleSimulator .
 **/
/*****************************************************************************/
class MoveBase
{
	public:
	//! Standard constructor (empty). Setting the current probability to Unity.
    MoveBase():probability(1.0){}


	/**
	 * @brief Update the acceptance probability e.g. multiply current probability with factor.
	 *
	 * @param factor Another probability of the move to update with current probability.
	 **/
	void multiplyProbability(double factor){probability*=factor;}

	/**
	 * @brief Returns the current acceptance probability
	 *
	 * @return current acceptance probability
	 **/
	double getProbability() const {return probability;}

	/**
	 * @brief Reset the current acceptance probability to 1.0
	 **/
	void resetProbability(){probability=1.0;}

private:
	//! Current probability of the move to be accepted by a Metropolis-criterion. See FeatureBoltzmann.
	double probability;

};


#endif
