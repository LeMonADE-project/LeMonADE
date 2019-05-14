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

#ifndef LEMONADE_IO_ABSTRACTWRITE_H
#define LEMONADE_IO_ABSTRACTWRITE_H

#include <iostream>

/*****************************************************************************/
/**
 * @file
 * @brief Definition of classes SuperAbstractWrite and AbstractWrite
 * */
/*****************************************************************************/


/*****************************************************************************/
/**
 * @class SuperAbstractWrite
 * @brief Abstract base class for writing commands
 *
 * @var headerOnly
 * @brief if true, the command should only execute once at the beginning
 * */
/*****************************************************************************/
class SuperAbstractWrite
{
public:
  SuperAbstractWrite():headerOnly(false){}
  virtual ~SuperAbstractWrite(){}

  virtual void writeStream(std::ostream& strm) = 0;
  bool writeHeaderOnly(){return headerOnly;}
protected:
  bool headerOnly;
};

/*****************************************************************************/
/**
 * @class AbstractWrite
 * @brief Extends SuperAbstractWrite by a reference to the source of information to be written
 * */
/*****************************************************************************/
template < class Source >
class AbstractWrite : public SuperAbstractWrite
{
public:
  AbstractWrite(const Source& src):src(src){}
  virtual ~AbstractWrite(){}

  virtual void writeStream(std::ostream& strm) {};


  const Source& getSource() const {return src;}
protected:
  void setHeaderOnly(bool state){headerOnly=state;}
private:
  const Source& src;


};

#endif

