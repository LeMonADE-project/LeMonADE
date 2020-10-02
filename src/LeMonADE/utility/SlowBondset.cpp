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

#include <LeMonADE/utility/SlowBondset.h>

/***********************************************************/
/**
 * @file
 * @brief Definition of methods of class SlowBondset
 */
/***********************************************************/

using namespace Lemonade;

SlowBondset::SlowBondset():bondsetLookup(0),lookupOffset(0){}

/***********************************************************/
/**
 * @brief Copy constructor for SlowBondset
 *
 * @param copyBondSet
 */
SlowBondset::SlowBondset(const SlowBondset& copyBondSet):
  FastBondset(), lookupOffset(copyBondSet.lookupOffset)
	    {
#ifdef DEBUG
		// check call of copy constructor
	        std::cout << "Bondset copy constructor" << std::endl;
#endif /*DEBUG*/
		//copy bondset
		this->BondVectors=copyBondSet.BondVectors;

	        // perform the lookup update
		  // set manually field pointer to 0 to avoid memory leak
		this->bondsetLookup=NULL;
		  // perform the reset
		this->updateLookupTable();
	    }

SlowBondset::~SlowBondset(){

	if(bondsetLookup!=0){
	    int32_t lookupSize=2*(lookupOffset)+1;
	    for(int32_t i=0;i<lookupSize;++i){
	      for(int32_t j=0;j<lookupSize;++j){
		delete[] bondsetLookup[i][j];
	      }
	      delete[] bondsetLookup[i];
	    }
	    delete[] bondsetLookup;
	  }
}


/***********************************************************/
/**
 * @details Before adding the bond vector given as argument to the set,
 * the function checks if both vector and id are still available,
 * i.e. not used for a different bondvector in the set.
 *
 * @throw <std::runtime_error> if the \a _bondVector or \a _identifier are already used.
 *
 * @param _bondVector Bond-Vector to add to the set.
 * @param _identifier Corresponding AscII-identifier for the bond.
 */
void SlowBondset::addBond(VectorInt3 _bondVector, int32_t _identifier)
{

  bool vectorExists=false;
  bool idExists=false;

  //check if bondvector already exists in the set
  for(iterator it=BondVectors.begin();it!=BondVectors.end();++it){

    if(it->second==_bondVector){
      vectorExists=(it->second==_bondVector);
      break;
    }
    if(it->first==_identifier){
      idExists=(it->first==_identifier);
      break;
    }
  }

  //add the bond, if it does not exist already, throw exception otherwise
  if(!(idExists||vectorExists)){
    std::cout<<" accepted bond with id "<<_identifier<<std::endl;
    BondVectors[_identifier]=_bondVector;
  }
  else{

    std::stringstream errormessage;

    errormessage<<"SlowBondset::addBond(VectorInt3 bondVector):\n"
	<<"bond vector with components "
	<<_bondVector.getX()<<" "<<_bondVector.getY()<<" "<<_bondVector.getZ()
	<<"\n or identifier "<<_identifier<<" already existent";

    throw std::runtime_error(errormessage.str());
  }

  //lookup table is now out of sync
  lookupSynchronized=false;
}

/**
 * @details Before adding the bond vector specified by the arguments to the set,
 * the function checks if both vector and identifier are still available,
 * i.e. not used for a different bondvector in the set.
 *
 * @throw <std::runtime_error> if the x,y,z-components exceed -4 or 4
 *
 * @throw <std::runtime_error> if the \a _bondVector or \a _identifier are already used.
 *
 * @param x The Cartesian x-component of the bond-vector.
 * @param y The Cartesian y-component of the bond-vector.
 * @param z The Cartesian z-component of the bond-vector.
 * @param identifier Corresponding AscII-identifier for the bond.
 */
void SlowBondset::addBond(int32_t x, int32_t y, int32_t z, int32_t identifier)
{
  VectorInt3 temp(x,y,z);
  addBond(temp, identifier);
  //lookup table is now out of sync
  lookupSynchronized=false;
}


/**
 * @details If the set of bond-vectors and the lookup table are
 * out of sync, the lookup is newly constructed. If the set of bond-vectors
 * is empty, a dummy value bondsetLookup[0][0][0]=false is produced.
 **/
void SlowBondset::updateLookupTable()
{
  //execute only if out of sync
  if(!lookupSynchronized){
    resetLookupTable();

    std::cout<<"setting up bondset lookup\n";

    //if no bonds are set, the lookup holds only a dummy value at 0,0,0
    if(size()==0){
      bondsetLookup=new bool**;
      bondsetLookup[0]=new bool*;
      bondsetLookup[0][0]=new bool;
      bondsetLookup[0][0][0]=false;
      lookupSynchronized=true;
      return;
    }
    else{
      //get max allowed bond-vector coordinate (absolute value)
      iterator it;
      for(it=begin();it!=end();++it){
	//positive values
	if( it->second.getX() > lookupOffset ) lookupOffset=it->second.getX();
	if( it->second.getY() > lookupOffset ) lookupOffset=it->second.getY();
	if( it->second.getZ() > lookupOffset ) lookupOffset=it->second.getZ();
	//negative values
	if( (-1)*(it->second.getX()) > lookupOffset ) lookupOffset=(-1)*(it->second.getX());
	if( (-1)*(it->second.getY()) > lookupOffset ) lookupOffset=(-1)*(it->second.getY());
	if( (-1)*(it->second.getZ()) > lookupOffset ) lookupOffset=(-1)*(it->second.getZ());
      }
      //now add one to the offset to make the lookup larger than the allowed vectors by at least 1 in every direction
      lookupOffset++;
      //determine the size of the lookup array in every direction (+1 vor x=y=z=0)
      int32_t lookupSize=2*(lookupOffset)+1;
      //now allocate memory for bondsetLookup and set all vectors to false
      bondsetLookup=new bool**[lookupSize];

      for(int32_t i=0;i<lookupSize;++i){
	bondsetLookup[i]=new bool*[lookupSize];

	for(int32_t j=0;j<lookupSize;++j){
	  bondsetLookup[i][j]=new bool[lookupSize];

	  for(int32_t k=0;k<lookupSize;++k){
	    bondsetLookup[i][j][k]=false;
	  }
	}

      }

      //now set the existing vectors true
      for(it=begin();it!=end();++it){
	bondsetLookup[it->second.getX()+lookupOffset]
		      [it->second.getY()+lookupOffset]
		      [it->second.getZ()+lookupOffset]=true;
      }
    }
    //lookup table is now synchronized with set of bondvectors
    lookupSynchronized=true;
  }
}

// Deletes the look-up table
void SlowBondset::resetLookupTable()
{
  if(bondsetLookup!=0){
    int32_t lookupSize=2*(lookupOffset)+1;
    for(int32_t i=0;i<lookupSize;++i){
      for(int32_t j=0;j<lookupSize;++j){
	delete[] bondsetLookup[i][j];
      }
      delete[] bondsetLookup[i];
    }
    delete[] bondsetLookup;
    bondsetLookup=0;
  }
  lookupOffset=0;
  lookupSynchronized=false;
}


/**
 * @throw <std::runtime_error> if look-up table is not synchronized with the system.
 *
 * @param bondVector Reference to VectorInt3 as bond-vector to check.
 * @return True if bond-vector is allowed, false otherwise.
 */
bool SlowBondset::isValid(const VectorInt3& bondVector ) const
{
  //lookup must not be out of sync
  if(!lookupSynchronized) throw std::runtime_error("SlowBondset::isValid(VectorInt3&)  Lookup not synchronized");
  //shift coordinates by lookupOffset
  //using uint here to make sure negative values turn into high positive ones,
  //such that the condition below works
  uint32_t x=bondVector.getX()+lookupOffset;
  uint32_t y=bondVector.getY()+lookupOffset;
  uint32_t z=bondVector.getZ()+lookupOffset;

  //if bondvectors are within coordinate range of lookup, see if they are valid
  //lookup offset cannot be negative, so conversion is save
  if(x<=2*uint32_t(lookupOffset) && y<=2*uint32_t(lookupOffset) && z<=2*uint32_t(lookupOffset))
    return bondsetLookup[x][y][z];
  //otherwise always return false
  else return false;
}

/**
 * @throw <std::runtime_error> if look-up table is not synchronized with the system.
 *
 * @param bondVector Reference to VectorInt3 as bond-vector to check.
 * @return True if bond-vector is allowed, false otherwise.
 */
bool SlowBondset::isValidStrongCheck(const VectorInt3& bondVector ) const
{
  return isValid(bondVector);
}
