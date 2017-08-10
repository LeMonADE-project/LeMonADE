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

#ifndef LEMONADE_FEATURE_FEATUREBOXREAD_H
#define LEMONADE_FEATURE_FEATUREBOXREAD_H

#include <stdexcept>

#include <LeMonADE/io/AbstractRead.h>

/**
 * @file
 * @brief Reading routines for bfm Reads !box_x,!box_y,!box_z,!periodic_x,!periodic_y,!periodic_z
 *
 */

//forward declarations
template <class IngredientsType> class ReadBoxX;
template <class IngredientsType> class ReadBoxY;
template <class IngredientsType> class ReadBoxZ;
template <class IngredientsType> class ReadPeriodicX;
template <class IngredientsType> class ReadPeriodicY;
template <class IngredientsType> class ReadPeriodicZ;



/*******************************************************************
 * Definitions of Read classes for the different bfm Reads
 *******************************************************************/


/*****************************************************************/
/**
 * @class ReadBoxX
 *
 * @brief Handles BFM-File-Reads \b !box_x
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class ReadBoxX:public AbstractRead
{
public:
  ReadBoxX(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};

/***********************************************************************/
/**
 * @class ReadBoxY
 *
 * @brief Handles BFM-File-Reads \b !box_y
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class ReadBoxY:public AbstractRead
{
public:
  ReadBoxY(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};

/***********************************************************************/
/**
 * @class ReadBoxZ
 *
 * @brief Handles BFM-File-Reads \b !box_z
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class ReadBoxZ:public AbstractRead
{
public:
  ReadBoxZ(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};

/***********************************************************************/
/**
 * @class ReadPeriodicX
 *
 * @brief Handles BFM-File-Reads \b !periodic_x
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class ReadPeriodicX:public AbstractRead
{
public:
  ReadPeriodicX(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};

/***********************************************************************/
/**
 * @class ReadPeriodicY
 *
 * @brief Handles BFM-File-Reads \b !periodic_y
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class ReadPeriodicY:public AbstractRead
{
public:
  ReadPeriodicY(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};

/***********************************************************************/
/**
 * @class ReadPeriodicZ
 *
 * @brief Handles BFM-File-Reads \b !periodic_z
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class ReadPeriodicZ:public AbstractRead
{
public:
  ReadPeriodicZ(IngredientsType& system):bfmSystem(system){};
  void execute();
private:
  IngredientsType& bfmSystem;
};


/******************************************************************
 * execute()-member function definitions of Read objects
 *****************************************************************/

/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !box_x
 *
 * @throw <std::runtime_error> if box_x could not be read.
 **/
template <class IngredientsType>
void ReadBoxX <IngredientsType>::execute()
{
	std::cout<<"reading box x dimension...";
  int x;
  if(*source>>x){
    bfmSystem.setBoxX(x);
    std::cout<<x<<std::endl;
  }
  else
    throw std::runtime_error("ReadBoxX::execute()\nCould not read box x dimension");
}


/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !box_y
 *
 * @throw <std::runtime_error> if box_y could not be read.
 **/
template <class IngredientsType>
void ReadBoxY <IngredientsType>::execute()
{
	std::cout<<"reading box y dimension...";
  int y;
  if(*source>>y){
    bfmSystem.setBoxY(y);
    std::cout<<y<<std::endl;
  }
  else
    throw std::runtime_error("ReadBoxY::execute()\nCould not read box y dimension");
}


/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !box_z
 *
 * @throw <std::runtime_error> if box_z could not be read.
 **/
template <class IngredientsType>
void ReadBoxZ <IngredientsType>::execute()
{
	std::cout<<"reading box z dimension...";
  int z;
  if(*source>>z){
    bfmSystem.setBoxZ(z);
    std::cout<<z<<std::endl;
  }
  else
    throw std::runtime_error("ReadBoxZ::execute()\nCould not read box z dimension");
}


/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !periodic_x
 *
 * @throw <std::runtime_error> if !periodic_x could not be read.
 **/
template <class IngredientsType>
void ReadPeriodicX <IngredientsType>::execute()
{
	std::cout<<"reading box x periodicity...";
  bool temp;
  if(*source>>temp){
    bfmSystem.setPeriodicX(temp);
    std::cout<<temp<<std::endl;
  }
  else
    throw std::runtime_error("ReadPeriodicX::execute()\nCould not read box x periodicity");
}


/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !periodic_y
 *
 * @throw <std::runtime_error> if !periodic_y could not be read.
 **/
template <class IngredientsType>
void ReadPeriodicY <IngredientsType>::execute()
{
	std::cout<<"reading box y periodicity...";
  bool temp;
  if(*source>>temp){
    bfmSystem.setPeriodicY(temp);
    std::cout<<temp<<std::endl;
  }
  else
    throw std::runtime_error("ReadPeriodicY::execute()\nCould not read box y periodicity");
}


/***********************************************************************/
/**
 * @brief Executes the reading routine to extract \b !periodic_z
 *
 * @throw <std::runtime_error> if !periodic_z could not be read.
 **/
template <class IngredientsType>
void ReadPeriodicZ <IngredientsType>::execute()
{
	std::cout<<"reading box z periodicity...";
  bool temp;
  if(*source>>temp){
    bfmSystem.setPeriodicZ(temp);
    std::cout<<temp<<std::endl;
  }
  else
    throw std::runtime_error("ReadPeriodicZ::execute()\nCould not read box z periodicity");
}

#endif
