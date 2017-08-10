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

#ifndef LEMONADE_CORE_GENERATECONTEXTTYPE_H
#define LEMONADE_CORE_GENERATECONTEXTTYPE_H

/**
 * @file
 *
 * @brief Classes used to generate Ingredients base class from features.
 **/

#include "extern/loki/Typelist.h"
#include "extern/loki/HierarchyGenerators.h"

#include <LeMonADE/core/TypelistExtensions.h>
#include <LeMonADE/core/FeatureHolder.h>
#include <LeMonADE/feature/Feature.h>

/**
 * @class EmptyModelFeature
 * @brief Used as dead end of the feature hierarchy by GenerateContextType (defined below).
 **/
struct EmptyModelFeature : public Feature {

};


/**
 * @class FullExpand
 * @brief Fully expands the Feature list (template parameter) with all (sub-)dependencies.
 *
 **/
template <class FeatureList> struct FullExpand;

//! Specialization for the case of no Features defined
template <> struct FullExpand< ::Loki::NullType >{typedef ::Loki::NullType Result;};

/**
 * @brief Implementation for the normal case with a typelist of features
 */
template < class Head, class Tail >
struct FullExpand< ::Loki::Typelist< Head, Tail > >
{
  //put the fully expanded forwards required features of the acutal feature (Head)
  //in front of head
 /**
 * @typedef Result
 * @brief Loki typelist containing all features and their dependencies (including possible duplicates).
 *
 * @details The typedef is the typelist given as template parameter recursively extended by all
 * dependencies and sub-dependencies of the features in this list. The class is used by
 * InsertFeatureRequests. Duplicate features are possible in the result and are
 * deleted later by InsertFeatureRequests.
 **/
  typedef typename ::Loki::TL::AppendIfNotNull
  <
    //insert the actual feature's backwards required features and append thei
    //fully expanded rest of the typelist
    typename ::Loki::TL::Append
    <
      //typelist consisting of the actual feature (Head) followed by the fully expanded list
      //of the backwards required features of Head
      ::Loki::Typelist<Head,typename FullExpand<typename Head::required_features_back>::Result >
      //fully expanded rest of the typelist
      ,typename FullExpand<Tail>::Result

    >::Result
    //fully expanded forwards required features of Head
    ,typename FullExpand<typename Head::required_features_front>::Result
  >::Result Result;
};



/**
 * @class InsertBackRequests
 * @brief Inserts all backwards dependencies in a feature list (not recursively)
 *
 **/
//forward declaration
template< class TList> struct InsertBackRequests;

//! Specialization for the case of no features defined
template <> struct InsertBackRequests< ::Loki::NullType>{typedef ::Loki::NullType Result;};

//! Implementation for the normal case with a typelist of features
template< class Head, class Tail>
struct InsertBackRequests< ::Loki::Typelist <Head,Tail> >
{
	/**
	 * @typedef Result
	 * @brief Resulting extended typelist of features
	 *
	 * @details The typelist given as template parameter is extended by all features
	 * defined as required_features_back by any feature in the original list. In
	 * contrast to FullExpand this is not done recursively here. Possible duplicates
	 * in the result are removed by InsertFeatureRequests.
	 **/
  typedef typename ::Loki::TL::Append
  <
    ::Loki::Typelist <Head ,typename Head::required_features_back>
    ,typename InsertBackRequests<Tail>::Result
  >::Result Result;
};




/**
 * @class InsertFeatureRequests
 * @brief Inserts all forward and backward dependencies into a typelist of Features.
 *
 **/
template<class TList> struct InsertFeatureRequests
{
  //get rid of duplicates from front to back after applying InsertBackRequests
  typedef typename ::Loki::TL::NoDuplicatesReverse
  <
    //reinsert the backwards required features because the NoDuplicates may have
    //deleted backwards dependenci es
  /**
  * @typedef Result
  * @brief Resulting typelist including all dependencies and no duplicates.
  *
  * @details inserts all requested features in the required order (in most cases...)
  * procedure: OrigList->FullExpand->NoDuplicates->InsertBackRequests->NoDuplicatesReverse
  **/
    typename InsertBackRequests
    <
      //remove duplicates back to front from fully expanded list of features
      typename ::Loki::TL::NoDuplicates< typename FullExpand < TList >::Result >::Result
    >::Result
  >::Result Result;
};


/**
 * @class GenerateContextType
 * @brief Generates a linear hierarchy type from the typelist of features given as template parameter
 **/
template < class TList > struct  GenerateContextType
{
	typedef ::Loki::GenLinearHierarchy < TList , FeatureHolder, EmptyModelFeature >  Result;
};


#endif /*LEMONADE_CORE_GENERATECONTEXTTYPE_H*/
