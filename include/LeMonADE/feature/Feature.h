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

#ifndef LEMONADE_FEATURE_FEATURE_H
#define LEMONADE_FEATURE_FEATURE_H

#include <string>
#include <ostream>

#include "extern/loki/Typelist.h"
#include "extern/loki/TypelistMacros.h"
#include "extern/loki/Functor.h"
#include "extern/loki/HierarchyGenerators.h"

#include <LeMonADE/core/TypelistExtensions.h>
#include <LeMonADE/updater/moves/MoveBase.h>

/**
 * @file
 *
 * @class Feature
 *
 * @brief Base Feature class, which every special Feature is derived from
 *
 **/
class Feature
{
public:

	Feature() {};

	virtual ~Feature() {};
  /**
   * @typedef required_features_front;
   * @brief Other special Features on which the actual one depends on. As default this list is empty.
   *
   * @details This behavior is necessary if the specialized Feature needs a certain Feature.
   * For example: FeatureExcludedVolume (excluded volume) need FeatureLattice (provides a lattice).
   */
  typedef ::Loki::NullType required_features_front;

  /**
   * @typedef required_features_back
   * @brief Other special Features, whose existence is implied by the actual Feature.  As default this list is empty.
   *
   * @details This behavior is necessary if the specialized Feature implies a certain Feature.
   * For example: thermal interaction implies Boltzmann criterion FeatureBoltzmann.
   * If a Feature provides an energy difference in the Monte-Carlo-Step an Metropolis-criterion
   * has to be applied. After collecting all energy differences the move has to be evaluate e.g. by FeatureBoltzmann.
   *
   */
  typedef ::Loki::NullType required_features_back;

  /**
   * @typedef monomer_extensions
   * @brief A list of monomer properties, which are directly related to the actual Feature. As default this list is empty.
   *
   * @details The monomer is then decorated with an additional interface
   * (e.g. charge in case of electrostatic interaction feature, an AttributeTag in FeatureAttributes)
   *
   */
  typedef ::Loki::NullType monomer_extensions;

  //! Export the relevant functionality for reading bfm-files to the responsible reader object
  template < class FileRead  > void exportRead ( FileRead & ){}

  //! Export the relevant functionality for writing bfm-files to the responsible writer object
  template < class FileWrite > void exportWrite( FileWrite& ) const {}


  /**
   * @brief Check for all Move. Does Nothing - Return True for all implementations.
   *
   * @details This dummy function is implemented for generality and inheritance.
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system
   * @param [in] move General move (maybe MoveLocalSc or MoveLocalBcc).
   * @return true Always!
   */
  template < class IngredientsType> bool checkMove( const IngredientsType&, const MoveBase& ) const{ return true; }

  /**
   * @brief Apply all Move. Does Nothing.
   *
   * @details This function applies all moves.
   * It does nothing and is implemented for generality and inheritance.
   *
   * @param [in] ingredients A reference to the IngredientsType - mainly the system
   * @param [in] move General move (maybe MoveLocalSc or MoveLocalBcc).
   */
  template < class IngredientsType> void applyMove( IngredientsType& ingredients, const MoveBase& move ) { }

  /**
   * @brief Synchronizes the Features and establishing consistency in the system.
   *
   * @details This behavior is necessary regarding the actual Feature, to make the whole
   * system of Ingredients consistent, and to check also for consistency.
   * (e.g. all bond-vectors between monomers have to be checked by the FeatureBondset)
   * (e.g. lattice has to be filled by excluded volume feature FeatureExcludedVolume)
   * It does nothing and is implemented for generality and inheritance.
   *
   * @param ingredients A reference to the IngredientsType - mainly the system.
   */
  template < class IngredientsType > void synchronize(IngredientsType& ingredients) {};

  /**
   * @brief Overloaded function to stream all metadata to an output stream.
   *
   * @details It does nothing and is implemented for generality and inheritance.
   * But in concrete implementation it should give additional metadata comments.
   *
   * @param stream output stream
   */
  virtual void printMetaData(std::ostream& stream) const{}

};


#endif /* LEMONADE_FEATURE_FEATURE_H */
