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

#ifndef LEMONADE_FEATURE_FEATUREBOXWRITE_H
#define LEMONADE_FEATURE_FEATUREBOXWRITE_H

#include <LeMonADE/io/AbstractWrite.h>

/***********************************************************************/
/**
 * @file
 * @brief writing routines of FeatureBox
 * */
/***********************************************************************/


/***********************************************************************/
/**
 * @class WriteBoxX
 *
 * @brief Handles BFM-File-Write \b !box_x
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBoxX:public AbstractWrite<IngredientsType>
{
public:
  //! Only writes \b !box_x into the header of the bfm-file.
  WriteBoxX(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !box_x.
  void writeStream(std::ostream& strm){strm<<"!box_x="<<(this->getSource()).getBoxX()<<"\n\n";}
private:

};

/***********************************************************************/
/**
 * @class WriteBoxY
 *
 * @brief Handles BFM-File-Write \b !box_y
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBoxY:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !box_y into the header of the bfm-file.
  WriteBoxY(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !box_y.
  void writeStream(std::ostream& strm){strm<<"!box_y="<<(this->getSource()).getBoxY()<<"\n\n";}
private:

};

/***********************************************************************/
/**
 * @class WriteBoxZ
 *
 * @brief Handles BFM-File-Write \b !box_z
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBoxZ:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !box_z into the header of the bfm-file.
  WriteBoxZ(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !box_z.
  void writeStream(std::ostream& strm){strm<<"!box_z="<<(this->getSource()).getBoxZ()<<"\n\n";}
private:

};

/***********************************************************************/
/**
 * @class WritePeriodicX
 *
 * @brief Handles BFM-File-Write \b !periodic_x
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WritePeriodicX:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !periodic_x into the header of the bfm-file.
  WritePeriodicX(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !periodic_x.
  void writeStream(std::ostream& strm){strm<<"!periodic_x="<<(this->getSource()).isPeriodicX()<<"\n\n";}
private:

};

/***********************************************************************/
/**
 * @class WritePeriodicY
 *
 * @brief Handles BFM-File-Write \b !periodic_y
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WritePeriodicY:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !periodic_y into the header of the bfm-file.
  WritePeriodicY(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !periodic_y.
  void writeStream(std::ostream& strm){strm<<"!periodic_y="<<(this->getSource()).isPeriodicY()<<"\n\n";}
private:

};

/***********************************************************************/
/**
 * @class WritePeriodicZ
 *
 * @brief Handles BFM-File-Write \b !periodic_y
 *
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WritePeriodicZ:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b !periodic_z into the header of the bfm-file.
  WritePeriodicZ(const IngredientsType& src):AbstractWrite<IngredientsType>(src){this->setHeaderOnly(true);}

  //! Executes the routine to write \b !periodic_z.
  void writeStream(std::ostream& strm){strm<<"!periodic_z="<<(this->getSource()).isPeriodicZ()<<"\n\n";}
private:

};

#endif
