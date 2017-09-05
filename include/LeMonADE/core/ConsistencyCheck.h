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

#ifndef LEMONADE_CORE_CONSISTENCYCHECK_H
#define LEMONADE_CORE_CONSISTENCYCHECK_H

#include "extern/loki/Typelist.h"

#include <LeMonADE/core/TypelistExtensions.h>
/****************************************************************************/
/**
 * @file
 *
 * @class ConsistencyCheck
 *
 * @brief Check a list of features for consistency regarding the required features (dependencies)
 *
 * @details Checks (from front to back) whether a feature in the list was required
 * to be further in front by some previous feature. For every feature in the list it
 * is checked whether this particular feature is contained in the list ForbiddenFeatures, i.e.
 * was required to be further in front.
 * The number of such inconsistencies are collected and the value of MY_ERRORSTATE is set accordingly.
 *
 * @tparam FeatureList
 * @brief List of features to be checked
 *
 * @tparam ForbiddenFeatures
 * @brief List of features that are required to be further in front in the list (this should be resolved automatically).
 * */
/****************************************************************************/

// forward declaration
template< class FeatureList, class ForbiddenFeatures> class ConsistencyCheck;


/**
 * @brief Specialization for the case when the end of the FeatureList is reached.
 *
 * @details In this case the first template parameter is NullType, because that's how typelists
 * are terminated
 */
template<class ForbiddenFeatures>
class ConsistencyCheck< ::Loki::NullType,ForbiddenFeatures>
{
private:
  enum {MY_VALUE= 0};
public:
  enum {MY_ERRORSTATE= MY_VALUE};

};



/**
 * @brief Specialization for the beginning of the typelist.
 *
 * @details At the beginning there are no features in the list of forbidden features (second parameter is NullType).
 *
 */
template< class Head, class Tail>
class ConsistencyCheck< ::Loki::Typelist< Head,Tail >, ::Loki::NullType>
{
private:
  typedef typename Head::required_features_front UpdatedForbiddenFeatures;
  enum {MY_VALUE=0};
public:
  enum {MY_ERRORSTATE=ConsistencyCheck<Tail,UpdatedForbiddenFeatures>::MY_ERRORSTATE};
};

//

/**
 * @brief Standard case at any point in a typelist.
 *
 * @details For every type in the list the required_features_front are appended
 * to the forbidden features and it is checked if the current feature is part of that list.
 * the number of such inconsistencies for this feature are counted by
 * \a MY_VALUE and the total number of inconsistencies of this and the following features
 * in the list are counted by \a MY_ERRORSTATE.
 */
template <class Head, class Tail, class ForbiddenFeatures>
class ConsistencyCheck< ::Loki::Typelist<Head,Tail> , ForbiddenFeatures >
{
private:

  typedef typename ::Loki::TL::NoDuplicates
  <
    typename ::Loki::TL::AppendIfNotNull <ForbiddenFeatures,typename Head::required_features_front>::Result
  >::Result UpdatedForbiddenFeatures;

  enum {MY_VALUE=::Loki::TL::IndexOf<ForbiddenFeatures,Head>::value==-1 ? 0 : 1};
public:

  enum {MY_ERRORSTATE= MY_VALUE + ConsistencyCheck<Tail,UpdatedForbiddenFeatures>::MY_ERRORSTATE};

};

#endif
