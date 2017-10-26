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

#ifndef LEMONADE_UPDATER_ABSTRACTUPDATER_H
#define LEMONADE_UPDATER_ABSTRACTUPDATER_H

/**
 * @file
 *
 * @class AbstractUpdater
 *
 * @brief Generic base class for all Updater. Provides the pure \a virtual interfaces
 * initialize(), execute(), and cleanup(). All Updater have to implement these functions.
 *
 * @details Since this class simply serves as a common base, the functions here don't do
 * 			anything particular. It's only provides interfaces for the TaskManager.
 *
 **/
class AbstractUpdater
{
public:
	//! Standard constructor (empty)
  AbstractUpdater(){}

  //! Standard destructor (empty)
  virtual ~AbstractUpdater(){}

  /**
   * @brief This function is called in \a every \a cycle in the TaskManager.
   *
   * @details It´s a pure virtual function for inheritance.
   * Use this function for repeating tasks (e.g. !mcs-Reads)
   *
   * @return Depends on their specific implementation.
   **/
  virtual bool execute() = 0;

  /**
   * @brief This function is called \a once in the beginning of the TaskManager.
   *
   * @details It´s a pure virtual function for inheritance.
   * Use this function for initializing tasks (e.g. init SDL)
   *
   **/
  virtual void initialize() = 0;

  /**
   * @brief This function is called \a once in the end of the TaskManager.
   *
   * @details It´s a pure virtual function for inheritance.
   * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
   *
   **/
  virtual void cleanup() = 0;
};


#endif /* LEMONADE_UPDATER_ABSTRACTUPDATER_H */
