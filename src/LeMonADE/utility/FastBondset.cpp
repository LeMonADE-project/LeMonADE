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

#include <LeMonADE/utility/FastBondset.h>

/***********************************************************/
/**
 * @file
 * @brief Definition of methods of class FastBondset
 */
/***********************************************************/

using namespace Lemonade;


FastBondset::FastBondset():lookupSynchronized(false){}

FastBondset::~FastBondset(){}

/***********************************************************/
/**
 * @brief Copy constructor for FastBondset
 *
 * @param copyBondSet
 */
FastBondset::FastBondset(const FastBondset& copyBondSet):// Copy the map Bondvectors explicitly
  BondVectors(copyBondSet.BondVectors),
  lookupSynchronized(false)
	    {
#ifdef DEBUG
		// check call of copy constructor
	        std::cout << "Bondset copy constructor" << std::endl;
#endif /*DEBUG*/
	        // perform the lookup update
		this->updateLookupTable();
	    }


/**
 * @param x The Cartesian x-component of the bond-vector.
 * @param y The Cartesian y-component of the bond-vector.
 * @param z The Cartesian z-component of the bond-vector.
 *
 * @throw <std::runtime_error> if the bond-vector doesn't exist.
 *
 * @return Ascii code associated with this bond.
 */
int32_t FastBondset::getBondIdentifier(int32_t x, int32_t y, int32_t z) const
{

	VectorInt3 temp(x,y,z);

	//check the set of stored vectors for a vector with components (x,y,z)
	for(iterator it=BondVectors.begin();it!=BondVectors.end();++it)
	{
		if(it->second==temp) return it->first;
	}

	//if we are still there, the required bond does not exist in the set
	//throw an exception in this case
	std::stringstream errormessage;
	errormessage<<"FastBondset::getBondIdentifier: Bondvector "
		<<x<<" "<<y<<" "<<z<<" does not exist";
	throw std::runtime_error(errormessage.str());
}


/**
 * @param identifier Ascii code identifier associated with a bond.
 *
 * @throw <std::runtime_error> if the identifier doesn't exist.
 *
 * @return The bondvector associated with identifier
 */
const VectorInt3 FastBondset::getBondVector(int32_t identifier) const
{

	//check the set of stored vectors for the identifier
	iterator it=BondVectors.find(identifier);

	//check the set of stored vectors for a vector with components (x,y,z)
	if(it != BondVectors.end()) return (it->second);

	//if we are still there, the required bond does not exist in the set
	//throw an exception in this case
	std::stringstream errormessage;
	errormessage<<"FastBondset::getBondVector: Bondvector with identifier "
		<<identifier<<" does not exist";
	throw std::runtime_error(errormessage.str());
}



/**
 * @details Before adding the bond vector given as argument to the set,
 * the function checks if both vector and id are still available,
 * i.e. not used for a different bondvector in the set.
 *
 * @throw <std::runtime_error> if the x,y,z-components exceed -4 or 4
 *
 * @throw <std::runtime_error> if the \a _bondVector or \a _identifier are already used.
 *
 * @param _bondVector Bond-Vector to add to the set.
 * @param _identifier Corresponding AscII-identifier for the bond.
 */
void FastBondset::addBond(VectorInt3 _bondVector, int32_t _identifier)
{

	//if any of the components are larger than four, this bondset
	//cannnot handle the bond due to the method of checking
	//bonds for validity. If bonds with larger components are needed, one
	//has to use the slower, but more general class SlowBondset
	if(	_bondVector.getX()>4 ||
		_bondVector.getX()<-4||
		_bondVector.getY()>4 ||
		_bondVector.getY()<-4||
		_bondVector.getZ()>4 ||
		_bondVector.getZ()<-4)
	{
		std::stringstream errormessage;
		errormessage<<"FastBondset::addBond(VectorInt3 _bondVector, int32_t _identifier)\n"
				<<"with _bondVector "<<_bondVector<<" and identifier "<<_identifier<<std::endl
				<<"Cannot add a bond with component larger 4 or smaller -4, because\n"
				<<"of the method of checking for valid bonds in function FastBondset::isValid()\n"
				<<"If you have to use such bond vectors, use the class SlowBondset instead of\n"
				<<"class FastBondset. It is slower, but more generally usable.\n";
		throw std::runtime_error(errormessage.str());
	}

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
	if(!(idExists||vectorExists))
	{
		std::cout<<" accepted bond with id "<<_identifier<<std::endl;
		BondVectors[_identifier]=_bondVector;
	}
	else
	{

		std::stringstream errormessage;

		errormessage<<"FastBondset::addBond(VectorInt3 bondVector):\n"
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
void FastBondset::addBond(int32_t x, int32_t y, int32_t z, int32_t identifier)
{
	VectorInt3 temp(x,y,z);
	addBond(temp, identifier);
	//lookup table is now out of sync
	lookupSynchronized=false;
}


/**
 * @details If the set of bond-vectors and the look-up table are
 * out of sync, the look-up is newly constructed. If the set of bond-vectors
 * is empty, a dummy value bondsetLookup[0][0][0]=false is produced.
 */
void FastBondset::updateLookupTable()
{
	//execute only if out of sync
	if(!lookupSynchronized)
	{
		resetLookupTable();

		std::cout<<"setting up bondset lookup\n";

		//lookup table is now synchronized with set of bondvectors
		iterator it;
		for(it=BondVectors.begin();it!=BondVectors.end();++it)
		{
			bondsetLookup[bondVectorToIndex(it->second)]=true;
		}
		lookupSynchronized=true;
	}
}


void FastBondset::resetLookupTable()
{
	for(size_t n=0;n<512;n++) bondsetLookup[n]=false;
	lookupSynchronized=false;


}

void FastBondset::clear()
{
	BondVectors.clear();

	for(size_t n=0;n<512;n++) bondsetLookup[n]=false;
	lookupSynchronized=false;

}


// add the classic bccBFM bond-set to the map (50 vectors)
void FastBondset::addBccBFMclassicBondset()
{

	addBond(-2,0,0,65);
	addBond(0,-2,0,66);
	addBond(0,0,-2,67);
	addBond(2,0,0,68);
	addBond(0,2,0,69);
	addBond(0,0,2,70);
	addBond(-2,-2,0,71);
	addBond(-2,0,-2,72);
	addBond(0,-2,-2,73);
	addBond(-2,2,0,74);
	addBond(-2,0,2,75);
	addBond(0,-2,2,76);
	addBond(2,-2,0,77);
	addBond(2,0,-2,78);
	addBond(0,2,-2,79);
	addBond(2,2,0,80);
	addBond(2,0,2,81);
	addBond(0,2,2,82);
	addBond(-2,-2,-2,83);
	addBond(-2,-2,2,84);
	addBond(-2,2,-2,85);
	addBond(2,-2,-2,86);
	addBond(-2,2,2,87);
	addBond(2,-2,2,88);
	addBond(2,2,-2,89);
	addBond(2,2,2,90);
	addBond(-3,-1,-1,97);
	addBond(-3,-1,1,98);
	addBond(-3,1,-1,99);
	addBond(3,-1,-1,100);
	addBond(-3,1,1,101);
	addBond(3,-1,1,102);
	addBond(3,1,-1,103);
	addBond(3,1,1,104);
	addBond(-1,-3,-1,105);
	addBond(-1,-3,1,106);
	addBond(-1,3,-1,107);
	addBond(1,-3,-1,108);
	addBond(-1,3,1,109);
	addBond(1,-3,1,110);
	addBond(1,3,-1,111);
	addBond(1,3,1,112);
	addBond(-1,-1,-3,113);
	addBond(-1,-1,3,114);
	addBond(-1,1,-3,115);
	addBond(1,-1,-3,116);
	addBond(-1,1,3,117);
	addBond(1,-1,3,118);
	addBond(1,1,-3,119);
	addBond(1,1,3,120);
}


// add the classic scBFM bond-set to the map (108 vectors)
void FastBondset::addBFMclassicBondset()
{
	addBond(2,0,0,17);
	addBond(0,0,2,18);
	addBond(0,2,0,19);
	addBond(-2,0,0,20);
	addBond(0,0,-2,21);
	addBond(0,-2,0,22);

	addBond(2,1,0,23);
	addBond(1,0,2,24);
	addBond(0,2,1,25);
	addBond(-2,1,0,26);
	addBond(1,0,-2,27);
	addBond(0,-2,1,28);

	addBond(2,-1,0,29);
	addBond(-1,0,2,30);
	addBond(0,2,-1,31);
	addBond(-2,-1,0,32);
	addBond(-1,0,-2,33);
	addBond(0,-2,-1,34);

	addBond(1,2,0,35);
	addBond(2,0,1,36);
	addBond(0,1,2,37);
	addBond(-1,2,0,38);
	addBond(2,0,-1,39);
	addBond(0,-1,2,40);

	addBond(1,-2,0,41);
	addBond(-2,0,1,42);
	addBond(0,1,-2,43);
	addBond(-1,-2,0,44);
	addBond(-2,0,-1,45);
	addBond(0,-1,-2,46);

	addBond(2,1,1,47);
	addBond(1,1,2,48);
	addBond(1,2,1,49);
	addBond(-2,1,1,50);
	addBond(1,1,-2,51);
	addBond(1,-2,1,52);

	addBond(2,-1,1,53);
	addBond(-1,1,2,54);
	addBond(1,2,-1,55);
	addBond(-2,-1,1,56);
	addBond(-1,1,-2,57);
	addBond(1,-2,-1,58);

	addBond(2,1,-1,59);
	addBond(1,-1,2,60);
	addBond(-1,2,1,61);
	addBond(-2,1,-1,62);
	addBond(1,-1,-2,63);
	addBond(-1,-2,1,64);

	addBond(2,-1,-1,65);
	addBond(-1,-1,2,66);
	addBond(-1,2,-1,67);
	addBond(-2,-1,-1,68);
	addBond(-1,-1,-2,69);
	addBond(-1,-2,-1,70);

	addBond(3,0,0,71);
	addBond(0,0,3,72);
	addBond(0,3,0,73);
	addBond(-3,0,0,74);
	addBond(0,0,-3,75);
	addBond(0,-3,0,76);

	addBond(2,2,1,77);
	addBond(2,1,2,78);
	addBond(1,2,2,79);
	addBond(-2,2,1,80);
	addBond(2,1,-2,81);
	addBond(1,-2,2,82);
	addBond(2,-2,1,83);
	addBond(-2,1,2,84);
	addBond(1,2,-2,85);
	addBond(-2,-2,1,86);
	addBond(-2,1,-2,87);
	addBond(1,-2,-2,88);

	addBond(2,2,-1,89);
	addBond(2,-1,2,90);
	addBond(-1,2,2,91);
	addBond(-2,2,-1,92);
	addBond(2,-1,-2,93);
	addBond(-1,-2,2,94);
	addBond(2,-2,-1,95);
	addBond(-2,-1,2,96);
	addBond(-1,2,-2,97);
	addBond(-2,-2,-1,98);
	addBond(-2,-1,-2,99);
	addBond(-1,-2,-2,100);

	addBond(3,1,0,101);
	addBond(1,0,3,102);
	addBond(0,3,1,103);
	addBond(-3,1,0,104);
	addBond(1,0,-3,105);
	addBond(0,-3,1,106);
	addBond(3,-1,0,107);
	addBond(-1,0,3,108);
	addBond(0,3,-1,109);
	addBond(-3,-1,0,110);
	addBond(-1,0,-3,111);
	addBond(0,-3,-1,112);

	addBond(1,3,0,113);
	addBond(3,0,1,114);
	addBond(0,1,3,115);
	addBond(-1,3,0,116);
	addBond(3,0,-1,117);
	addBond(0,-1,3,118);
	addBond(1,-3,0,119);
	addBond(-3,0,1,120);
	addBond(0,1,-3,121);
	addBond(-1,-3,0,122);
	addBond(-3,0,-1,123);
	addBond(0,-1,-3,124);
}

