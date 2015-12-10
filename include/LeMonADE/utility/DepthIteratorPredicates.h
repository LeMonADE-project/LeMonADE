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

#ifndef LEMONADE_UTILITY_DEPTHITERATORPREDICATES_H
#define LEMONADE_UTILITY_DEPTHITERATORPREDICATES_H

#include <set>

#include <LeMonADE/utility/Vector3D.h>
/**
 * @file
 *
 * @brief Some predicates (functors) that can be used with the function fill_connected_groups
 **/

/**
 * @struct belongsToLinearStrand
 *
 * @brief Functor providing the information if Vertex is no cross-linking point
 *  (functionality 0 < f <= 2)
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorBelongToLinearStrand
 **/
struct belongsToLinearStrand {
	template<class MoleculesType>
	bool operator()(const MoleculesType& m, int i)
			{
		return ((m.getNumLinks(i) > 0) && (m.getNumLinks(i) < 3));
	}
};

/**
 * @struct alwaysTrue
 *
 * @brief Functor returning always \a True.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorAlwaysTrue
 **/
struct alwaysTrue
{
	template<class MoleculesType>
	bool operator()(const MoleculesType& m, int i)
			{
		return true;
	}
};

/**
 * @struct hasBonds
 *
 * @brief Functor returning \a True if monomer has at least one bond.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorAlwaysTrue
 **/
struct hasBonds
{
	template<class MoleculesType>
	bool operator()(const MoleculesType& m, int i)
			{
		return m.getNumLinks(i) >0;
	}
};

/**
 * @class hasType
 *
 * @brief Functor providing the information the Vertex has an attribute tag.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorHasType
 **/
class hasType
{
public:
	hasType(int type){types.insert(type);}
	hasType(int type1,int type2){types.insert(type1); types.insert(type2);}
	hasType(std::vector<int> typeVector){for(size_t n=0;n<typeVector.size();n++) types.insert(typeVector[n]);}
	
	template<class MoleculesType>
	bool operator()(const MoleculesType& m, int i)
	{
		return types.count(m[i].getAttributeTag());
	}
private:
	std::set<int> types;
};

/**
 * @class notOfType
 *
 * @brief Functor providing the information the Vertex does not have an attribute tag.
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo Rename FunctorNotOfType
 **/
template<int AttributeType>
class notOfType
{
public:
	
	template<class MoleculesType>
	bool operator()(const MoleculesType& m, int i)
	{
		return (m[i].getAttributeTag()==AttributeType?false:true);
	}

};

#endif
