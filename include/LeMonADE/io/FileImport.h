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

#ifndef LEMONADE_IO_FILEIMPORT_H
#define LEMONADE_IO_FILEIMPORT_H

/**
 *@file
 *
 *@todo proper reaction to unknown Read in function executeRead
 */

#include <string>
#include <fstream>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <stdint.h>
#include <stdexcept>

#include <LeMonADE/Version.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/Parser.h>



/*******************************************************************
 * definition of class template FileImport
 ******************************************************************/

/*********************************************************/
/**
 * @class FileImport
 *
 * @brief Manages the import of data from .bfm files
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 */
template <class IngredientsType>
class FileImport{

public:

  FileImport(const std::string& sourcefile,IngredientsType& dataStorage);
  virtual ~FileImport();


  //! Interface for adding new read Reads provided by Feature (e.g. !box_x)
  void registerRead(std::string, AbstractRead*);

  //! Interface for replacing old (standard) Read with new read Reads. Example: FeatureJumps (!mcs)
  void replaceRead(std::string, AbstractRead*);

  AbstractRead* modifyRead(std::string);

  void initialize()
  {
	  //read the header of the bfm file
	  readHeader();
	  scanFile();

	  //bfmData.synchronize(bfmData);
  }

  //! Parse and process bfm-file up to first conformation (incl. read-in first !mcs)
  void readHeader();

  //! Read, parse, and process file until eof or next !mcs is reached (or file stream fails otherwise)
  bool read();

  //! Jumps to the given !mcs in the file and reads it (unsafely forward-winding).
  bool gotoMcs(uint64_t mcs);

  //! Jumps to the given frame in the file and reads it (unsafely forward/reverse-winding).
  bool gotoFrame(uint32_t frame);


  //! Jumps to the to the first conformation in the file and reads it (unsafely backward-winding).
  bool gotoStart();

  //! Jumps to the last conformation in the file and reads it (unsafely forward-winding).
  bool gotoEnd();

  //! Reading file to the last conformation and reads it (safely forward-winding).
  bool gotoEndSave();

  //! Returns the smallest time (in mcs) in the file.
  uint64_t getMinAge();

  //! Returns the highest time (in mcs) in the file.
  uint64_t getMaxAge();

  //! Returns the number of frames (occurrence of !mcs) in the file.
  uint32_t getNumFrames(){return framePositionInFile.size();}



  //! Returns the recent frames (occurrence of !mcs) in the file.
  uint32_t getRecentFrameCount()
  {
	  ///todo is there any other possibility to get this information?
	  /*std::map <uint32_t, uint64_t>::const_iterator it;
	  for(it=framePositionInFile.begin();it!=framePositionInFile.end();++it)
	  {
		  if(it->second==bfmData.getMolecules().getAge()){
			  break;
		  }
	  }

	  return it->first;
	  */
	  return std::distance(mcsPositionInFile.begin(),mcsPositionInFile.find(bfmData.getMolecules().getAge()) )+1;
  }

  /**
   * @brief Returns the filename used in the FileImporter
   *
   * @return Filename+Suffix to read from.
   */
  std::string getFilename() const{return filename;}

  /**
   * @brief Set the filename (with suffix) used in the FileImporter.
   * @param filename File+Suffix to read from.
   */
  void setFilename(std::string filename){this->filename = filename;}

  //! Close the file stream.
  void close(){file.close();};

  //! Get pointer to data container
  IngredientsType& getDestination(){return bfmData;}

private:

  //! Storage for data that are read-in from file (mostly Ingredients).
  IngredientsType& bfmData;

  //file handling
  //! Name of bfm-file with file suffix *.bfm
  std::string filename;

  //! File stream associated with the input file (FileImport::file)
  std::ifstream file;

  //! Parser that finds and returns Read strings from input file
  Parser parser;

  //Read handling
  //! Map of Read-strings (e.g. !box_x) associated Read objects
  std::map <std::string, AbstractRead*> Reads;
  //save the unknown commands found in the file in
  //this set. This is used to make sure information about unknown reads is
  //not written out multiple times when printMetaData is used.
  //! Set of Read-strings of unknown commands found in the file.
  std::set <std::string> unknownReads;

  //! Execute read of Read associated with string
  void executeRead(const std::string&);

  // these are used for jumping to mcs-positions in the file
  //! The time (in mcs) of the first conformation.
  uint64_t firstMcs;

  //! Maps mcs to the position in the file right before the mcs.
  std::map<uint64_t,std::streampos> mcsPositionInFile;

  //! Maps frame number to mcs.
  std::map<uint32_t,uint64_t> framePositionInFile;

  //saves the complete first conformation. this is useful for
  //jumping back to the first conformation
  //! Holing the first conformation in the file. Useful for rewinding.
  typename IngredientsType::molecules_type firstConformation;

  //! Scans the complete file for the positions of the !mcs commands
  void scanFile();

};

/***********************************************************************
 * definitions of the methods of FileImport
 **********************************************************************/

/******************************************************
 *FileImport::FileImport(const string&,IngredientsType&)
 *****************************************************/
/**
 * @brief Constructor for basic file handling that reads header (incl 1.st conformation) and synchronizes all data.
 *
 * @details It initialize all internal values and passes all information to the corresponding classes.
 * It checks if the given file is ready for read-in operation. It register all known commands from
 * the Feature, reading the header of the file, and scans for all commands. At the end it synchronizes
 * itself with the Ingredients.
 *
 * @throw <std::runtime_error> No File Access.
 * @param sourcefile Name of the file to read in
 * @param dataStorage Class holding all information of the system (mainly Ingredients )
 *
 * @todo Did I understand that correctly that the first !mcs is read not the last in the file?
 */
template <class IngredientsType>
FileImport<IngredientsType>::FileImport(const std::string& sourcefile,IngredientsType& dataStorage)
  :bfmData(dataStorage),filename(sourcefile),parser(file),firstMcs(0)
{
  //open the source file
  file.open(sourcefile.c_str(),std::ios_base::in|std::ios_base::binary);
  if(file.fail()) throw std::runtime_error(std::string("error opening input file ")+sourcefile+std::string("\n"));

  //the features' read-Reads are registered here!!
  bfmData.exportRead(*this);

  //make the input file stream known to Read objects
  for(typename std::map <std::string,AbstractRead*>::iterator ib=Reads.begin();ib!=Reads.end();++ib){
    ib->second->setInputStream(&file);
  }

  dataStorage.setName(sourcefile);

  //read the header of the bfm file
  //readHeader();
  //scanFile();

  //bfmData.synchronize(bfmData);
}


/******************************************************
 *FileImport::~FileImport()
 *****************************************************/
template <class IngredientsType>
FileImport<IngredientsType>::~FileImport()
{
	std::map <std::string, AbstractRead*>::iterator it;
	for(it=Reads.begin();it!=Reads.end();++it){
		delete it->second;
		it->second=0;
	}
	Reads.clear();
}

/******************************************************************************
 *void FileImport::registerRead(string ReadString, AbstractRead* ReadObject)
 ******************************************************************************/
/**
 *
 * @param ReadString bfm-keyword/command/user-command (e.g. !box_x, #!fixed_monomers)
 * @param ReadObject pointer to object providing the Read's functionality
 */
template <class IngredientsType>
void FileImport<IngredientsType>::registerRead(std::string ReadString, AbstractRead* ReadObject)
{
  std::cout<<"registered bfm-Read "<<ReadString<<std::endl;
  if( Reads.insert(std::pair<std::string,AbstractRead*> (ReadString, ReadObject)).second==false)
  {
    std::stringstream errormessage;
    errormessage<<"Could not register Read with command string "<<ReadString
		<<" because the command string is already in use.";
    throw std::runtime_error(errormessage.str());
  }
}


/******************************************************************************
 *void FileImport::replaceRead(string key, AbstractRead* newRead)
 ******************************************************************************/
/**
 *
 * @param key bfm-keyword/command/user-command (e.g. !box_x, #!fixed_monomers)
 * @param newRead replace the read object with keyword key by the new one pointing to newRead
 */
template <class IngredientsType>
void FileImport<IngredientsType>::replaceRead(std::string key,AbstractRead* newRead)
{
	//find the old read. if not present, throw exception. otherwise replace
	std::map <std::string, AbstractRead*>::iterator it (this->Reads.find(key));
	if ( it == Reads.end() )
	{
		std::ostringstream strm; strm << "Read command \"" << key << "\" not found for modification.";
		throw std::runtime_error(strm.str());
	}
	else
	{
		delete it->second;
		it->second=newRead;
		std::cerr << "WARNING: REPLACED READ COMMAND \"" << key << "\".\n";
	}
}

///@todo KOMMENTIEREN
/**
 * @todo I don't get this? Thats the purpose of this? Thats the difference between replaceRead? It´s not used anywhere...
 *
 */
template <class IngredientsType> AbstractRead*
FileImport<IngredientsType>::modifyRead(std::string key)
{
	std::map <std::string, AbstractRead*>::iterator it (this->Reads.find(key));
	if ( it == Reads.end() )
	{
		std::ostringstream strm; strm << "Read command \"" << key << "\" not found for modification.";
		throw std::runtime_error(strm.str());
	}
	else
	{
		std::cout << "returning read command for modification \"" << key << "\".\n";
	}
	return it->second;
}

//Parse and process bfm-file up to first conformation (incl. read-in first !mcs)
template <class IngredientsType>
void FileImport<IngredientsType>::readHeader()
{
	double version;
	bool versionInformationPresent=false;

	std::string Read;
	std::streampos beforeFirstMcs=file.tellg();

	bool first_mcs_found=false;

	//read file up to the first mcs command (including the first mcs)
	while(!file.fail() && (Read != "endoffile") && (first_mcs_found == false))
	{
		Read=parser.findRead();


		if (Read=="!mcs")
		{
			if (!first_mcs_found) first_mcs_found = true;

		}
		else if(Read=="#!version")
		{
			file>>version;
			versionInformationPresent=true;

		}
		else
			beforeFirstMcs=file.tellg();

		executeRead(Read);

	}

	firstMcs=bfmData.getMolecules().getAge();
	//now reset the file to the position befor the first mcs command
	file.seekg(beforeFirstMcs);
	firstConformation=bfmData.getMolecules();

	if(!versionInformationPresent || version!=::LEMONADE_VERSION)
	{
		std::cerr<<"\nWARNING: NO VERSION INFORMATION PRESENT IN FILE, OR VERSION NUMBER DIFFERS FROM THE PROGRAM USED!\n";
	}
	std::cout << "readHeader : done " << std::endl;
}

/******************************************************************************
 *void FileImport::read()
 ******************************************************************************/
// read and process file until eof or next !mcs is reached (or file stream fails otherwise)
/**
 * @details If the header is not read-in before it also parses the header and the first !mcs.
 *
 * @return True if another !mcs is found. False if no !mcs is found after the last !mcs, anymore.
 */
template <class IngredientsType>
bool FileImport<IngredientsType>::read()
{

	std::string Read;

	//read and process file until eof or !mcs is reached (or file stream fails otherwise)
	bool MCSFound = false;

	while(!file.fail() && Read != "endoffile" )
	{
		Read=parser.findRead();
		if ( Read == "!mcs")
			MCSFound = true;

		executeRead(Read);

		if ( MCSFound ) break;
	}


	return MCSFound;
}

/******************************************************************************
 *void FileImport::executeRead(const string& ReadString)
 ******************************************************************************/
// Execute read of Read associated with string
/**
 * @param ReadString Keyword/command/user-command which should be executed.
 */
template <class IngredientsType>
void FileImport<IngredientsType>::executeRead(const std::string& ReadString)
{
  std::map<std::string,AbstractRead*>::iterator it;
  //try to find ReadString in the map of registered Reads
  it=Reads.find(ReadString);

  //if not found, react appropriately
  if (it==Reads.end()){
    //std::cout<<"unknown Read "<<ReadString<<std::endl;
    if(unknownReads.count(ReadString)==0)
    {
      bfmData.addComment(ReadString + "\n");
      unknownReads.insert(ReadString);
    }
    return;
  }
  //if found, call the appropriate function
  else
  {
	//std::cout << "found Read " << ReadString << std::endl;
    it->second->execute();
  }
}


/**
 * @details Scans the file for !mcs-commands and saves the positions of these commands in
 * the map mcsPositionInFile. The positions saved there are right after the
 * previous !mcs. It does \b NOT parse the !mcs nor it sets any system informations.
 *
 * @throw <runtime_error> if IO-error (different from eof) occurs.
 *
 * @todo Rename to scanFileForMCS()!
 **/
template <class IngredientsType>
void FileImport<IngredientsType>::scanFile()
{
	std::cout<<"scanning file...this may take some seconds for large files...";
	std::cout.flush();
	//save this position and return to it after the operation
	std::streampos startingPosition=file.tellg();

	std::string read;
	std::streampos mcsPosition;
	uint64_t mcs;

	//first go to the beginning of the file and find the position of
	//the first mcs command. this position is not saved, because it
	//is not necessarily clear, which commands preceeding the first
	// !mcs area part of this !mcs (e.g. solvent). Therefore, the
	//complete first conformation is saved in readHeader().
	file.clear();
	mcsPositionInFile.clear();
	framePositionInFile.clear();
	file.seekg(0,std::ios::beg);

	//the information of the first frame/conformation is read-in by readHeader()
	while(!file.fail() && (read != "endoffile"))
	{
		mcsPosition=file.tellg();

		read=parser.findRead();

		if (read=="!mcs")
		{
			file>>mcs;
			mcsPositionInFile.insert(std::make_pair(mcs,mcsPosition));
			framePositionInFile.insert(std::make_pair(mcsPositionInFile.size(),mcs));
			break;
		}
	}

	//now save all the following positions of the !mcs commands

	while(!file.fail() && read != "endoffile" )
	{

		mcsPosition=file.tellg();

		read=parser.findRead();
		while(read!="!mcs" && !file.fail() && read!="endoffile" )
			{
			read=parser.findRead();

			}

		if(!(read=="endoffile")) file>>mcs;

		if(file.fail() && !(read=="endoffile"))
		{
			std::stringstream errormessage;
            errormessage<<"FileImport::scanFile(): error scanning mcs positions after mcs "
                          <<mcsPositionInFile.rbegin()->first<<"\n";
			throw std::runtime_error(errormessage.str());
		}
		//std::cout<<"insterting mcs "<<mcs<<" at filepointer pos "<<mcsPosition<<std::endl;
		mcsPositionInFile.insert(std::make_pair(mcs,mcsPosition));
		framePositionInFile.insert(std::make_pair(mcsPositionInFile.size(),mcs));
	}

	//go back to the position the file was at before this function was called
	file.clear();
	file.seekg(startingPosition);

	std::cout<<"done\n";

}

//jumps to the mcs given as argument and reads the conformation
/**
 * @details IMPORTANT: TOPOLOGY MIGHT NOT BE CORRECT IF IT IS CHANGING SOMEWHERE IN THE
 * MIDDLE OF THE FILE! It parses the time smaller or equal to the required one.
 * It returns if there is another !mcs after the required one.
 *
 * @param mcs The time in MCS in the file for forward-winding
 * @return True if another !mcs is there. False, otherwise.
 */
template <class IngredientsType>
bool FileImport<IngredientsType>::gotoMcs(uint64_t mcs)
{
	//reset the eof-bit if its set
	file.clear();

  bool retVal;
  std::streampos position;
  std::map<uint64_t,std::streampos>::iterator it;
  //check if smaller than min or larger than max

  //if the required mcs is smaller than the first in the map, just
  //update the system with the first conformation
  if(mcs< mcsPositionInFile.begin()->first)
  {
    (bfmData.modifyMolecules())=firstConformation;
    file.seekg(mcsPositionInFile.begin()->second);
    retVal=true;
  }
  //similarily, if it is larger than the largest mcs contained in the file,
  //go to the last conformation
  else if(mcs>mcsPositionInFile.rbegin()->first)
  {
    it=mcsPositionInFile.end();
    --it;
    position=it->second;
    file.seekg(position);
    retVal=read();
  }
  //else got to the conformation with mcs smaller or equal to the
  //one required
  else
  {
    it=mcsPositionInFile.lower_bound(mcs);
    position=it->second;
    file.seekg(position);
    retVal=read();
  }

  return retVal;
}

/**
 * @details IMPORTANT: TOPOLOGY MIGHT NOT BE CORRECT IF IT IS CHANGING SOMEWHERE IN THE
 * MIDDLE OF THE FILE! It parses the time smaller or equal to the required one.
 * It returns if there is another !mcs after the required one.
 *
 * @param frame The frame in the file for forward/rewinding-winding
 * @return True if another !mcs is there. False, otherwise.
 */
template <class IngredientsType>
bool FileImport<IngredientsType>::gotoFrame(uint32_t frame)
{
	return gotoMcs(framePositionInFile.find(frame)->second);
}

//jumps to the last comformation in the file
/**
 * @details IMPORTANT: TOPOLOGY MIGHT NOT BE CORRECT IF IT IS CHANGING SOMEWHERE IN THE
 * MIDDLE OF THE FILE! It returns if there is another !mcs after last (should be false).
 *
 * @return True if another !mcs is there. False, otherwise.
 *
 * @todo Return value seems not to be necessary, do we?
 */
template <class IngredientsType>
bool FileImport<IngredientsType>::gotoEnd()
{

	std::streampos position;

  position=mcsPositionInFile.rbegin()->second;
  file.clear();

  file.seekg(position);
  //now read until next mcs is processed
  bool retVal=read();

  return retVal;
}

//goes through the file step by step until last conformation is reached
//this can be much slower than gotoEnd(), but it garantees that the
//topology is correct, in case it changes somewhere in the middle of the file
/**
 * @details IMPORTANT: IN CONTRAST TO gotoEnd() THE TOPOLOGY SHOULD BE CORRECT,
 * BECAUSE THE FILE IS READ STEP BY STEP AND EVERY COMMAND SHOULD BE PARSED.
 *
 * @return True if the eof is found.
 */
template<class IngredientsType>
bool FileImport<IngredientsType>::gotoEndSave()
{
	gotoStart();

	bool reachedEnd;
	do
	{
		reachedEnd=!(read());

	}while(reachedEnd==false);

	return reachedEnd;
}

//goes back to the first conformation. the file pointer is set accordingly in
//gotoMcs(0)
/**
 * @details IMPORTANT: TOPOLOGY MIGHT NOT BE CORRECT IF IT IS CHANGING SOMEWHERE IN THE
 * MIDDLE OF THE FILE! It parses the the first conformation in the file.
 * It returns if there is another !mcs after the required one.
 *
 * @todo Implement an safety backward-winding, or is this not necessary?
 * @return True if another !mcs is there. False, otherwise.
 */
template<class IngredientsType>
bool FileImport<IngredientsType>::gotoStart()
{
	bool retVal=gotoMcs(0);
	return retVal;
}

//returns the smallest mcs in the file
/**
 *
 * @return smallest time (in !mcs) in the file
 */
template<class IngredientsType>
uint64_t FileImport<IngredientsType>::getMinAge()
{
    return firstMcs;
}

//returns the highest mcs in the file
/**
 * @throw <std::runtime_error> if no !mcs is in the file
 *
 * @return highest time (in !mcs) in the file
 */
template<class IngredientsType>
uint64_t FileImport<IngredientsType>::getMaxAge()
{
    if(!mcsPositionInFile.empty())
        return mcsPositionInFile.rbegin()->first;
    else
        throw std::runtime_error("FileImport::getMaxAge(): no Monte-Carlo steps in file.\n");
}


#endif
