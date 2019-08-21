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

#ifndef FEATURE_BENDING_POTENTIAL_READ_WRITE_H
#define FEATURE_BENDING_POTENTIAL_READ_WRITE_H

/**
 * @file
 * @date 2019/08/11
 * @author Ankush Checkervarty
 * @brief Def. and impl. of class templates ReadBendingPotential and WriteBendingPotential
**/

#include<iostream>
#include<LeMonADE/io/AbstractRead.h>
#include<LeMonADE/io/AbstractWrite.h>

/**
 * @class ReadBendingPotential
 * @brief Handles BFM-file read command !bending_potential
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template < class IngredientsType>
class ReadBendingPotential: public ReadToDestination<IngredientsType>
{
public:
    ReadBendingPotential(IngredientsType& i):ReadToDestination<IngredientsType>(i){
//         i.modifyBondset().addBFMclassicBondset();
    }
    virtual ~ReadBendingPotential(){}
    virtual void execute();
};

/**
 * @class WriteBendingPotential
 * @brief Handles BFM-file write command !bending_potential
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteBendingPotential:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the potential
  //is written only once at the beginning of the output file.
    WriteBendingPotential(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteBendingPotential(){}

    virtual void writeStream(std::ostream& strm);
};

/////////////MEMBER IMPLEMENTATIONS ////////////////////////////////////////////

/**
 * @brief Executes the reading routine to extract \b !bending_potential.
 *
 * @throw <std::runtime_error> fail to read monomer types or bending potential strength.
 **/
template<class IngredientsType>
void ReadBendingPotential<IngredientsType>::execute()
{
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();
    
    
    int32_t type;
    float bpStrength;

    //set stream to throw exception on fail
    file.exceptions(file.exceptions() | std::ifstream::failbit);

    try
      {
	file>>type;
	file>>bpStrength;
      }
    catch(std::ifstream::failure e)
      {
	std::stringstream errormessage;
	errormessage<<"ReadBendingPotential::execute().\n";
	errormessage<<"Could not read interaction from file\n";
	errormessage<<"Previous error: "<<e.what()<<std::endl;
	throw std::runtime_error(errormessage.str());
      }
      
    //now save the interaction tuple just read from the file
    ingredients.setBendingPotential(ingredients, type,bpStrength);

    //unset exception on fail
    file.exceptions(file.exceptions() & (~std::ifstream::failbit));
}


/**
 * @brief Executes the routine to write \b !bending_potential.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteBendingPotential<IngredientsType>::writeStream(std::ostream& stream)
{
  int32_t nSpecies=10; //number must fit into 8 bit (uint8_t based lattice)
  stream<<"## bending potential magnitude at type in kT (default 0.0kT)\n";

  for(int32_t type=1;type<=nSpecies;type++)
    {
  	  if(this->getSource().getBendingPotential(type)!=0.0)
            {
	      stream<<"!bending_potential "<<type<<" "<<this->getSource().getBendingPotential(type)<<"\n";
            }

        }
  stream<<"\n\n";

}




#endif // FEATURE_NN_INTERACTION_READ_WRITE_H
