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

#ifndef FEATURE_NN_INTERACTION_READ_WRITE_H
#define FEATURE_NN_INTERACTION_READ_WRITE_H

/**
 * @file
 * @date 2016/06/18
 * @author Hauke Rabbel
 * @brief Def. and impl. of class templates ReadNNInteraction and WriteNNInteraction
**/

#include<iostream>
#include<LeMonADE/io/AbstractRead.h>
#include<LeMonADE/io/AbstractWrite.h>

/**
 * @class ReadNNInteraction
 * @brief Handles BFM-file read command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template < class IngredientsType>
class ReadNNInteraction: public ReadToDestination<IngredientsType>
{
public:
    ReadNNInteraction(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadNNInteraction(){}
    virtual void execute();
};

/**
 * @class WriteNNInteraction
 * @brief Handles BFM-file write command !nn_interaction
 * @tparam IngredientsType Ingredients class storing all system information.
**/
template <class IngredientsType>
class WriteNNInteraction:public AbstractWrite<IngredientsType>
{
public:

  //constructor sets the headerOnly tag, such that the interaction
  //is written only once at the beginning of the output file.
    WriteNNInteraction(const IngredientsType& i)
        :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}

    virtual ~WriteNNInteraction(){}

    virtual void writeStream(std::ostream& strm);
};

/////////////MEMBER IMPLEMENTATIONS ////////////////////////////////////////////

/**
 * @brief Executes the reading routine to extract \b !nn_interaction.
 *
 * @throw <std::runtime_error> fail to read monomer types or interaction energy.
 **/
template<class IngredientsType>
void ReadNNInteraction<IngredientsType>::execute()
{
    IngredientsType& ingredients=this->getDestination();
    std::istream& file=this->getInputStream();

    int32_t typeA,typeB;
    double interactionConstant;

    //set stream to throw exception on fail
    file.exceptions(file.exceptions() | std::ifstream::failbit);

    try
      {
	file>>typeA;
	file>>typeB;
	file>>interactionConstant;
      }
    catch(std::ifstream::failure e)
      {
	std::stringstream errormessage;
	errormessage<<"ReadNNInteraction::execute().\n";
	errormessage<<"Could not read interaction from file\n";
	errormessage<<"Previous error: "<<e.what()<<std::endl;
	throw std::runtime_error(errormessage.str());
      }

    //now save the interaction tuple just read from the file
    ingredients.setNNInteraction(typeA,typeB,interactionConstant);

    //unset exception on fail
    file.exceptions(file.exceptions() & (~std::ifstream::failbit));
}


/**
 * @brief Executes the routine to write \b !nn_interaction.
 * @arg stream file stream to write into
 **/
template<class IngredientsType>
void WriteNNInteraction<IngredientsType>::writeStream(std::ostream& stream)
{
  int32_t nSpecies=255; //number must fit into 8 bit (uint8_t based lattice)
  stream<<"## nearest neighbor interactions between types in kT (default 0.0kT)\n";

  for(int32_t typeA=1;typeA<=nSpecies;typeA++)
    {
      for(int32_t typeB=1;typeB<=typeA;typeB++)
        {
	  if(this->getSource().getNNInteraction(typeA,typeB)!=0.0)
            {
	      stream<<"!nn_interaction "<<typeB<<" "<<typeA<<" "<<this->getSource().getNNInteraction(typeB,typeA)<<"\n";
            }

        }

    }
  stream<<"\n\n";

}




#endif // FEATURE_NN_INTERACTION_READ_WRITE_H
