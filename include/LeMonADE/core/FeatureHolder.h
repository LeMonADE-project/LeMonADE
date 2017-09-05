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

#ifndef LEMONADE_CORE_FEATUREHOLDER_H_
#define LEMONADE_CORE_FEATUREHOLDER_H_

#include <iostream>
#include <cxxabi.h>

#include "extern/loki/Typelist.h"
#include "extern/loki/TypelistMacros.h"
#include "extern/loki/Functor.h"
#include "extern/loki/HierarchyGenerators.h"

#include <LeMonADE/core/TypelistExtensions.h>
/**
 * @file
 *
 * @class FeatureHolder
 *
 * @brief Used by class GenerateContextType to assemble the simulated system from used Features.
 *
 * @details Besides this is being used by GenerateContextType to put together the simulated or
 * analyzed system from a list of features, this class also facilitates
 * communication and synchronization of all features. Functions common to all
 * features can thus be called for every used feature in the system (see
 * member functions)
 *
 * @tparam Feature
 * @brief Feature held by this FeatureHolder
 *
 * @tparam Base
 * @brief Linear hierarchy of the rest of the used features
 *
 **/
template < class Feature, class Base > class FeatureHolder: public Feature, public Base
{

public:

  	/**
	 * @brief Export the relevant functionality for reading bfm-files of all used Features.
	 *
	 * @param file File importer for the bfm-file
	 **/
  template < class FileRead  > void exportRead ( FileRead& file )
  {
	//register read functionalities of this feature
	Feature::exportRead(file);
	//register read functionalities of the rest of the features
	Base   ::exportRead(file);
  }

  /**
   * @brief Export the relevant functionality for writing bfm-files of all used Features.
   *
   * @param file File importer for the bfm-file
   **/
  template < class FileWrite > void exportWrite( FileWrite& file) const
  {
	//register write functionalities of this feature
	Feature::exportWrite(file);
	//register write functionalities of the rest of the features
	Base   ::exportWrite(file);
  }

  /**
   * @brief Overloaded function to stream all metadata of the Features to an output stream.
   *
   * @details It demangels the name the specified Feature by using @a abi::__cxa_demangle() and delegates
   * the stream to the features and base class.
   *
   * @param stream output stream
   */
  void printMetaData(std::ostream& stream) const
  {
	int status;
	char * demangled = abi::__cxa_demangle(typeid(Feature(*this)).name(),0,0,&status);
	stream << demangled << std::endl;
	//print meta data of this feature
	Feature::printMetaData(stream);
	//print meta data of the rest of the features
	Base   ::printMetaData(stream);
	//free memory
	free(demangled);
  }

  /**
   * @brief Check for all Features the Move.
   *
   * @details It delegates the check to all Features and Base.
   * Returns \a true if move is allowed or rejected (\a false ).
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system
   * @param [in] move General move
   * @return true if move is allowed or rejected (\a false ).
   */
  template < class IngredientsType, class MoveType > bool checkMove( const IngredientsType& ingedients, MoveType& move )
  {
    //check move compatibility with current feature and
    // the rest of the features.
    return Feature::checkMove(ingedients,move) && Base::checkMove(ingedients,move);
  }


  /**
   * @brief This function applies the given Move to the Ingredients
   *
   * @details It delegates the applyMove to all Features and Base.
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system
   * @param [in] move General move.
   */
  template < class IngredientsType, class MoveType > void applyMove( IngredientsType& ingedients, const MoveType& move )
  {
    //check move compatibility with current feature and
    // the rest of the features.
    Feature::applyMove(ingedients,move);
    Base::applyMove(ingedients,move);
  }


  /**
   * @brief Synchronizes all Features and establishing consistency in the system.
   *
   * @details It delegates the synchronization to all Features and Base
   * to make the Ingredients system consistent.
   *
   * @param ingredients A reference to the IngredientsType - mainly the system.
   */
  template < class IngredientsType > void synchronize(IngredientsType& ingredients)
  {
	Feature::synchronize(ingredients); Base::synchronize(ingredients);
  }

};

#endif /* LEMONADE_CORE_FEATUREHOLDER_H_ */
