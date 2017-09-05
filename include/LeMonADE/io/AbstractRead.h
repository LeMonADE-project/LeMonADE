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

#ifndef LEMONADE_IO_ABSTRACTREAD_H
#define LEMONADE_IO_ABSTRACTREAD_H

/***********************************************************************/
/**
 *@file
 *@brief definition of base class for reading operations from bfm-file
 */
/***********************************************************************/

#include <fstream>
#include <string>
#include <vector>
#include <sstream>

/***********************************************************************/
 /**
  * @class AbstractRead
  *
  * @brief Base class for reading Reads that are given by the Feature.
  *
  **/
class AbstractRead
{
public:

  //! Default constructor (empty)
  AbstractRead(){};

  //! Default destructor (empty)
  virtual ~AbstractRead(){};

  //! Pure virtual function for execution of read-in
  virtual void execute()=0;

  /**
   * @brief Sets the input stream by the given parameter \a stream.
   *
   * @param stream Specified input stream
   */
  void setInputStream(std::istream* stream){source=stream;};

protected:

  /**
   * @brief Return the used input stream.
   *
   * @return used input stream.
   */
  std::istream& getInputStream(){return *source;};

  //! Convenience function for detecting a command line
  bool detectRead(std::string& line) const;

  //! Convenience function for detecting a separator character
  bool findSeparator(std::istream& stream,char separator);

  //! Convenience function for splitting up a string into substrings
  std::vector<std::string> tokenize2Parameter(const std::string& str,
  		char delim1, char delim2);

  //! Pointer to the stream the Read is reading in
  std::istream* source;

};


/**
 * @class ReadToDestination
 * @brief Extends the base class AbstractRead by a pointer to be a container for information
 *
 * */
template < class Destination > class ReadToDestination : public AbstractRead
{
  public:
	/**
	 * @brief Default constructor. Copy the given destination.
	 * @param destination A given reference to store the data
	 */
  ReadToDestination(Destination& destination):dst(destination){}

  //! Default destructor (empty)
  virtual ~ReadToDestination(){};

  /**
   * @brief Returns a reference to structure storing the data from read-in
   *
   * @return Structure stroing the data
   */
  virtual Destination& getDestination() {return dst;}

  private:
  //! Pointer to the structure that stores the data read from file
  Destination& dst;
};

#endif
