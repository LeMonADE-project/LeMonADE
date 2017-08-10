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

#ifndef LEMONADE_CORE_GENERATEMONOMERTYPE_H
#define LEMONADE_CORE_GENERATEMONOMERTYPE_H

/**
 * @file
 * @brief class templates NextMonomerExtension, MonomerExtensionRequests and GenerateMonomerType
 **/

#include "extern/loki/Typelist.h"
#include "extern/loki/TypelistMacros.h"
#include "extern/loki/Functor.h"
#include "extern/loki/HierarchyGenerators.h"

#include <LeMonADE/core/TypelistExtensions.h>
#include <LeMonADE/core/FeatureHolder.h>

/**
 * @struct NextMonomerExtension
 *
 * @brief Finds the next non-NullType monomer extension in the list of features.
 *
 * @details Assigns a typedef Result to either the template parameter MonoFeature, or, if
 * MonoFeature is of ::Loki::NullType, assigns this typedef to the next non-NullType
 * monomer extensions in the remaining list of Features. Also assigns a typedef
 * Rest to the remaining list of features following the next non-NullType
 * monomer extension. This class is used by MonomerExtensionRequests to create a
 * typelist of necessary monomer extensions.
 *
 * @tparam MonoFeature
 * @brief Monomer_extension(s) of the current Feature in the list to be searched
 *
 * @tparam TList
 * @brief tail of the list of Features following the current one
 **/
template <class MonoFeature,class TList> struct NextMonomerExtension
{
  typedef MonoFeature Result;
  typedef TList Rest;
};

//! Specialization for the case of MonoFeature being of ::Loki::NullType
template<class Head, class Tail>
struct NextMonomerExtension< ::Loki::NullType, ::Loki::Typelist<Head,Tail> >
{
    typedef typename NextMonomerExtension<typename Head::monomer_extensions,Tail>::Result Result;
    typedef typename NextMonomerExtension<typename Head::monomer_extensions,Tail>::Rest Rest;
};

//! Specialization for the case of no requests at all
template<> struct NextMonomerExtension< ::Loki::NullType, ::Loki::NullType>
{
  typedef  ::Loki::NullType Result;
  typedef  ::Loki::NullType Rest;
};


/******************************************************************************/

/**
 * @struct MonomerExtensionRequests
 * @brief Creates a typelist of monomer extensions required by the features in TList.
 *
 * @tparam TList
 * @brief List of Features to be searched for monomer extensions.
 *
 * @typedef Result
 * @brief Resulting typelist of monomer extensions.
 **/
template <class TList> struct MonomerExtensionRequests;

//! Specialization for the case of no requests at all
template<> struct MonomerExtensionRequests < ::Loki::NullType >
{
  typedef ::Loki::NullType Result;
};

//! Normal case with a typelist<Head,Tail> to be searched for required extensions
template<class Head, class Tail>
struct MonomerExtensionRequests< ::Loki::Typelist<Head,Tail> >
{
  typedef typename ::Loki::TL::NoDuplicates
  <
      typename ::Loki::TL::AppendIfNotNull
      <
	typename MonomerExtensionRequests< typename NextMonomerExtension< typename Head::monomer_extensions,Tail>::Rest >::Result,
	typename NextMonomerExtension< typename Head::monomer_extensions,Tail>::Result
      > ::Result

  > ::Result Result;
};

/******************************************************************************/

/**
 * @struct AddFeaturesIfNotNull
 * @brief Helper class for GenerateMonomerType
 **/
template<class MonomerBase, class TList>
struct AddFeaturesIfNotNull{
  typedef ::Loki::GenLinearHierarchy < TList, FeatureHolder, MonomerBase> Result;
};

//! Specialization for the case of empty TList
template<class MonomerBase>
struct AddFeaturesIfNotNull<MonomerBase,::Loki::NullType>{
  typedef MonomerBase Result;
};

/******************************************************************************/

/**
 * @struct GenerateMonomerType
 * @brief Decorates MonomerBase with the features listed in the typelist TList, if TList is not empty. Otherwise
 * it leaves MonomerBase as it is
 **/
template<class MonomerBase, class TList>
struct GenerateMonomerType
{
  typedef typename AddFeaturesIfNotNull<MonomerBase,typename MonomerExtensionRequests<TList>::Result>::Result Result;

};


#endif  /*LEMONADE_CORE_GENERATEMONOMERTYPE_H*/
