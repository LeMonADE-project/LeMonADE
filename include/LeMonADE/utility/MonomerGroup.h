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

#ifndef LEMONADE_UTILITY_MONOMERGROUP_H
#define LEMONADE_UTILITY_MONOMERGROUP_H

#include <stdexcept>
#include <sstream>

/**
 * @brief Basic operation on MonomerGroup
 *
 * @deprecated
 *
 * @todo we should reconsider this approach for usability
 *
 * @todo delete?
 */
template < class MoleculesType > class MonomerGroup
{

  const MoleculesType* moleculesGroup;

  std::vector < int > indices;

public:

  typedef typename MoleculesType::vertex_type vertex_type;

  MonomerGroup( const MoleculesType& moleculesGroup_):moleculesGroup(&moleculesGroup_){}

  void push_back(int idx){indices.push_back(idx);}

  //! removes the n-th particle from the group. n is NOT the true index of the particle!
  void erase(int n);

  //! removes the monomer with index n from the group. here n is the true index of the particle in MoleculesType m
  void removeFromGroup(int n);

  const vertex_type& operator[] (int i) const { return (*moleculesGroup)[indices.at(i)];}

  void operator += (const MonomerGroup& rhs) {
    //loop over all elemtns of right hand side
    for(size_t i=0;i<rhs.indices.size();i++){
      std::vector<int>::iterator it;
      for(it=indices.begin();it!=indices.end();++it){
	//check if element is already there
	if(*it == rhs.indices.at(i)) break;
      }
      if(it==indices.end()){
	// if not, add element
	indices.push_back(rhs.indices.at(i));
      }
    }
  }

  MonomerGroup& operator = (const MonomerGroup& src) {
    this->moleculesGroup = src.moleculesGroup;
    this->indices = src.indices;
    return *this;
  }

  size_t size() const { return indices.size();}

  int trueIndex( int i ) const {return indices.at(i);}

  uint64_t getAge() const {return moleculesGroup->getAge();}

  MoleculesType copyGroup() const;

  void clear(){indices.clear();}

};


template<class MoleculesType>
MoleculesType MonomerGroup<MoleculesType>::copyGroup() const
{
	//make a temporary molecules object
	MoleculesType mol;

	//set age and size
	mol.resize(size());
	mol.setAge(getAge());

	//get a map of all indices (as in the original molecules) to the indices
	//in the group
	std::map<int,int> mapOrigToGroupIndex;
	for(size_t n=0;n<size();++n)
	{
		mapOrigToGroupIndex.insert( std::pair<int,int>(indices.at(n),n));
	}
	//copy particles and links
	//first loop over all particles in the group
	for(size_t n=0;n<size();++n)
	{
		//stores the index the n-th particle of the group has in
		//the original molecules
		int originalIndex=indices.at(n);

		//copy this monomer to the new molecules
		mol[n]=(*moleculesGroup)[originalIndex];

		//now copy the bonds and bond information
		//loop over all bond partners of the current monomer
		for(size_t l=0; l<moleculesGroup->getNumLinks(originalIndex); l++)
		{
			//index of the l-th neighbor in original molecules
			int origNeighborIndex=moleculesGroup->getNeighborIdx(originalIndex,l);

// 			//if that neighbor is also part of this group, copy the link
			if( mapOrigToGroupIndex.count(origNeighborIndex)>0 && mol.areConnected(n,mapOrigToGroupIndex.at(origNeighborIndex))==false )
			{
				mol.connect(n,mapOrigToGroupIndex.at(origNeighborIndex),moleculesGroup->getLinkInfo(originalIndex,origNeighborIndex));
			}
		}

	}

	//return the newly constructed molecules object
	return mol;

}

template<class MoleculesType>
void MonomerGroup<MoleculesType>::erase(int n)
{

	if(n<int(indices.size()) && n>=0)
		indices.erase(indices.begin()+n);
	else{

		std::stringstream errormessage;
		errormessage<<"MonomerGroup::erase(int n): n="<<n<<" index out of bounds";
		throw std::runtime_error(errormessage.str());
	}
}

template<class MoleculesType>
void MonomerGroup<MoleculesType>::removeFromGroup(int n)
{
	if(indices.size()==0){
		std::stringstream errormessage;
		errormessage<<"MonomerGroup::removeFromGroup(int n): n="<<n<<"...cannot delete,group is empty";
		throw std::runtime_error(errormessage.str());
	}

	std::vector<int>::iterator it;

	for(it=indices.begin();it!=indices.end();++it){
		if(*it == n) break;
	}

	if(it!=indices.end()) indices.erase(it);
	else{
		std::stringstream errormessage;
		errormessage<<"MonomerGroup::removeFromGroup(int n): n="<<n<<" is not part of the group";
		throw std::runtime_error(errormessage.str());
	}
}


#endif /* LEMONADE_UTILITY_MONOMERGROUP_H */
