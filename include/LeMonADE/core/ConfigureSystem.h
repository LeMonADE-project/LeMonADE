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

#ifndef LEMONADE_CORE_CONFIGURESYSTEM_H
#define LEMONADE_CORE_CONFIGURESYSTEM_H

#include "extern/loki/Typelist.h"

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/GenerateContextType.h>
#include <LeMonADE/core/GenerateMonomerType.h>
#include <LeMonADE/core/ConsistencyCheck.h>

/**
 * @file
 *
 * @class ConfigureSystem
 *
 * @brief Collection of typedefs for configuring Ingredients from a list of features.
 *
 * @tparam max_connectivity Maximum connectivity of the monomers. Default is 7.
 *
 * @tparam Edge Type of information that can be stored on a bond. Default is \a int.
 *
 * @typedef ConfigureSystem::feature_list
 * @brief Final list of features used to construct Ingredients (Loki typelist).
 * @details Here the list of features given as a template parameter is automatically
 * extended by all features that may be required by other features in the list.
 *
 * @typedef ConfigureSystem::context_type
 * @brief Base type for Ingredients
 * @details The type consists of a linear inheritance hierarchy of all features in
 * the typelist feature_list, connected by FeatureHolder.
 *
 * @typedef ConfigureSystem::monomer_type
 * @brief Basic type for monomers used in Molecules.
 * @details The MonomerBaseType given as the first template argument is extended by
 * all monomer extensions required by features in feature_list.
 *
 * @typedef ConfigureSystem::molecules_type
 * @brief Defines the graph used to store the monomers.
 * @details Monomers of the previously defined monomer_type are here stored in Molecules,
 * which provides interfaces and functionalities for handling the monomers.
 *
 * @enum MY_ERRORSTATE
 * @brief 0 if everything good, > 0 else
 * @details Counts the number of inconsistencies in feature_list. More precisely
 * this means the number of features that should be further towards the head of
 * the list as required by other features' dependencies.
 * */


template <class MonomerBaseType, class FeatureList, uint max_connectivity=7, class Edge=int> class ConfigureSystem
{
public:
  typedef typename InsertFeatureRequests < FeatureList >::Result feature_list;
  typedef typename GenerateContextType<feature_list>::Result context_type;
  typedef typename GenerateMonomerType<MonomerBaseType,feature_list>::Result monomer_type;
  typedef Molecules<monomer_type,max_connectivity,Edge> molecules_type;
  enum {MY_ERRORSTATE=ConsistencyCheck< feature_list,::Loki::NullType>::MY_ERRORSTATE};
private:
};

#endif
