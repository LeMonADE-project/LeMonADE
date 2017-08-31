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

#ifndef LEMONADE_UPDATER_UPDATERREADBFMFILE_H
#define LEMONADE_UPDATER_UPDATERREADBFMFILE_H

/*****************************************************************************/
/**
 * @file
 * @brief Definition of Updater UpdaterReadBfmFile
 * */
/*****************************************************************************/

#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/io/FileImport.h>

/*****************************************************************************/
/**
 * @class UpdaterReadBfmFile
 *
 * @brief Updater that reads configurations from bfm-file. On each execute call it reads the next configuration.
 **/
template <class IngredientsType>
class UpdaterReadBfmFile: public AbstractUpdater
{
public:
	/**
	 * @brief Standard (empty) constructor initializing the variables
	 *
	 * @param filename PathAndFilename+Suffix to be loaded
	 * @param ing A reference to the IngredientsType - mainly the system
	 * @param readType ENUM for specify read-in
	 */
  UpdaterReadBfmFile(std::string filename,IngredientsType& ing,int readType)
    :ingredients(ing),file(filename,ing),myReadType(readType){};


  /**
   * @enum BFM_READ_TYPE
   *
   * @brief Specifiers for bfm-file-handling what should be read-in.
   */
  enum BFM_READ_TYPE{
  	READ_LAST_CONFIG_SAVE=0, //!< Read last configuration in file (safe mode)
  	READ_LAST_CONFIG_FAST=1, //!< Read last configuration in file (unsafe mode)
  	READ_STEPWISE=2			 //!< Read first configuration in file then stepwise (safe mode)
  };

  /**
   * @brief This function is called \a once in the beginning of the TaskManager.
   *
   * @details This function initializes the UpdaterReadBFM with the given conformation in the file.
   *          After succeeding the read-in depending which ENUM you choose, the system is synchronized
   *          with the first conformation/frame (READ_STEPWISE) or the last conformation/frame
   *          (READ_LAST_CONFIG_SAVE or READ_LAST_CONFIG_FAST).
   *
   *
   * @throw <std::runtime_error> If unknown ENUM-type is used.
   *
   **/
  virtual void initialize()
  {

	  file.initialize();


		if(myReadType==READ_STEPWISE)
			file.read();
		else if(myReadType==READ_LAST_CONFIG_SAVE)
		{
			file.gotoEndSave();
			//retVal=false;
		}
		else if(myReadType==READ_LAST_CONFIG_FAST)
		{
			file.gotoEnd();
			//retVal=false;
		}
		else
			throw std::runtime_error("ReadBFMFile::init()...read type must have value READ_LAST_CONFIG_SAVE(0), READ_LAST_CONFIG_FAST(1), or READ_STEPWISE(2)\n");

		ingredients.synchronize();
  };


  /**
   * @brief Read-in of the next configuration. Synchronize the Ingredients.
   *
   * @details This function executes the UpdaterReadBFM read-in the next conformation in the file.
   *          The system is synchronized with the read-in conformation/frame (READ_STEPWISE).
   *          This does nothing if the myReadType is READ_LAST_CONFIG_SAVE or READ_LAST_CONFIG_FAST.
   *          If end-of-file is reached it returns false, otherwise true.
   *
   *
   * @throw <std::runtime_error> If unknown ENUM-type is used.
   *
   * @return True if another !mcs is found. False if no !mcs is found after the last !mcs, anymore (EOF).
   */
  virtual bool execute()
  {
	bool retVal = false;

	if(myReadType==READ_STEPWISE)
		retVal=file.read();

	/*
	else if(myReadType==READ_LAST_CONFIG_SAVE)
	{
		file.gotoEndSave();
		retVal=false;
	}
	else if(myReadType==READ_LAST_CONFIG_FAST)
	{
		file.gotoEnd();
		retVal=false;
	}
	else
	*/
	if( ((myReadType != READ_STEPWISE) && (myReadType != READ_LAST_CONFIG_SAVE) && (myReadType!=READ_LAST_CONFIG_FAST)))
		throw std::runtime_error("ReadBFMFile::execute()...read type must have value READ_LAST_CONFIG_SAVE(0), READ_LAST_CONFIG_FAST(1), or READ_STEPWISE(2)\n");


	ingredients.synchronize();
	return retVal;
  }


  /**
   * @brief This function is called \a once in the end of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
   *
   **/
  virtual void cleanup(){};

  /**
   * @brief Returns the smallest time (in mcs) in the file.
   *
   * @return smallest time (in !mcs) in the file
   **/
  uint64_t getMinAge(){return file.getMinAge();}


  /**
   * @brief Returns the highest time (in mcs) in the file.
   *
   * @throw <std::runtime_error> if no !mcs is in the file
   *
   * @return highest time (in !mcs) in the file
   */
  uint64_t getMaxAge(){return file.getMaxAge();}

  /**
   * @brief Returns the number of frames (occurrence of !mcs)in the file.
   *
   * @throw <std::runtime_error> if no !mcs is in the file
   *
   * @return All number of frames in the file
   */
  uint32_t getNumFrames(){return file.getNumFrames();}

    /**
     * @brief Returns the recent frames (occurrence of !mcs) in the file.
     *
     * @throw <std::runtime_error> if no !mcs is in the file
     *
     * @return Recent frame in the file
     **/
  uint32_t getRecentFrameCount(){return file.getRecentFrameCount();}


  //! Jumps to the last conformation in the file and reads it (unsafely forward-winding).
  void gotoEnd(){file.gotoEnd();}

  //! Jumps to the last conformation in the file and reads it (unsafely reverse-winding).
  void gotoStart(){file.gotoStart();}

  //! Jumps to frame and reads conformation in the file and reads it (unsafely verse-winding).
  void gotoFrame(uint32_t frame){file.gotoFrame(frame);}


  //! Close the file stream.
  void closeFile(){file.close();}

private:
	//! A reference to the IngredientsType - mainly the system
	IngredientsType& ingredients;

	//! FileImport to read-in the configurations.
	FileImport<IngredientsType> file;

	//! ENUM-type BFM_READ_TYPE specify the read-in
	int myReadType;

};

#endif
