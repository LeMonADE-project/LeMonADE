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
 * @class Move
 *
 * @brief Base class for any types of moves inside the simulation.
 *
 * @details All moves must be derived from this base class in order to be processed
 * correctly by the Feature. All moves should override the functions contained here.
 * Since this class simply serves as a common base, the functions here don't do 
 * anything particular.
 **/
/*****************************************************************************/
class MoveBase
{
	public:
	//! Standard constructor (empty). Setting the current probability to Unity.
    MoveBase():probability(1.0){}

	/**
	* @brief Checks for acceptance of this move by the argument.
	*
	* @details This function hands the move to the ingredients object and
	* returns true, if all Features accept the move (and false otherwise). For most
	* cases this function can simply be copy-pasted into every specialized move as it is.
	* 
	* @param ingredients system to be checked for compatibility with the move
	* @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
	*
	* @todo This can be deleted because it´s never called. All Called are delegated to the SpecializedMove
	**/
	template<class IngredientsType>
	bool check(IngredientsType& ingredients) const
	{
	return ingredients.checkMove(ingredients,*this);   
	}
  
	/**
	* @brief Applies the move to the system Ingredients given as argument. Since this class simply serves as a common base, the function don't do
	* anything particular. Please read on details...
	*
	* @details When implementing this function one has to do two things: FIRST
	* apply the move to the system by calling ingredients.applyMove(..), THEN one
	* can do something specific to the move, e.g. change particle positions. Note 
	* the order of the two actions!!!
	*
	* @param ingredients system to be changed by the move
	*
	* @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
	*
	* @todo This can be deleted because it´s never called. All Called are delegated to the SpecializedMove
	**/
	template <class IngredientsType> void apply(IngredientsType& ingredients)
	{
	
		//all moves should then be applied to the features like this:
		//ingredients.applyMove(ingredients,*this);
		
	
		//specialized moves may do some particular update here, e.g. change particle
		//positions in some particular manner. For example:
		//ingredients.modifyMolecules()[1].setX(0);
	}
	
	/**
	* @brief Initializes the move. Since this class simply serves as a common base, the function don't do
	* anything particular. Please read on details...
	*
	* @details This is where a new random move should be drawn in the implementation.
	*
	* @param ingredients system for which the move is initialized
	* 
	* @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
	*
	* @todo This can be deleted because it´s never called. All Called are delegated to the SpecializedMove
	**/
	template <class IngredientsType> void init(const IngredientsType& ingredients){};
	
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
