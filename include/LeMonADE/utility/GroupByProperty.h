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

#ifndef LEMONADE_UTILITY_GROUPBYPROPERTY_H
#define LEMONADE_UTILITY_GROUPBYPROPERTY_H
/**
 * @file
 * @brief definition of global function group_by_property
 * */

#include <map>
#include <string>

/**
 *
 * @brief Global function that groups vertices/monomers by a property
 *
 * @details The function combines indexes of monomers into groups defined by a certain property. Each group is of type GroupType, and the
 * groups themselves are stored in the GroupContainer groups. For example, GroupContainer could be a vector of lists of indices.
 * The functor findGroup provides an ()-operator which is used to separate the monomers into groups. Its return value is the key
 * that identifies a group. For example, if one wanted to group monomers with respect to their number of bonds, the operator would
 * simply return the number of bonds of a certain monomer.
 *
 * @param molecules reference to the container/graph, in which the monomers are stored
 * @param groups reference to the container that stores the groups to be found here
 * @param groupInit reference to an object that serves as base type for a group (argument group stores objects of this type)
 * @param findGroup functor that separates the monomers into groups using its ()-operator. It needs to provide a typedef property_type to specify the return type.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 **/
template < class Graph, class GroupContainer, class GroupType, class GroupFinder >
void group_by_property(Graph& molecules, GroupContainer& groups, const GroupType& groupInit, const GroupFinder& findGroup)
{

  typedef typename GroupFinder::property_type Property;
  //all values of the property in question are mapped to an index (size_t)
  //this is the index of the group as used in the GroupContainer groups.
  std::map<Property,size_t> detectedGroups;

  //loop over all monomers to find get the property and put them into the correct group
  for(size_t i=0;i<molecules.size();++i)
  {
    //if the group key, returned by findGroup(i), already exists in the map,
    //append the index of the current monomer to the respective group
    if(detectedGroups.find(findGroup(molecules,i))!=detectedGroups.end()){
      groups[detectedGroups.find(findGroup(molecules,i))->second].push_back(i);
    }
    //if the group key does not exist yet, open up a new group in groups
    //and append the index of the current monomer to that group
    else{
      detectedGroups.insert(std::make_pair(findGroup(molecules,i),groups.size()));
      groups.push_back(groupInit);
      groups.back().push_back(i);
    }

  }
}

/**
 * @class NumberOfNeighbors
 *
 * @brief Example for a functor for use with group_by_property. It sorts the monomers by their number of bonds
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorNumberOfNeighbors
 **/
struct NumberOfNeighbors
{
  typedef int property_type;
  template<class Graph> property_type operator()(Graph& molecules,int i) const {return molecules.getNumLinks(i);}
};

/**
 * @class ZCoordinate
 *
 * @brief Example for a functor for use with group_by_property. It sorts the monomers into groups by their z-coordinate in intervals of 10.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorZCoordinate
 **/
struct ZCoordinate
{
  typedef int property_type;
  template<class Graph> property_type operator()(Graph& molecules, int i) const {return 10* (int)(molecules[i].getZ()/10);}
};

#endif
