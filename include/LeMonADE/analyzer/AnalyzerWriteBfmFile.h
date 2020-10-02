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

#ifndef LEMONADE_ANALYZER_ANALYZERWRITEBFMFILE_H
#define LEMONADE_ANALYZER_ANALYZERWRITEBFMFILE_H

#include <set>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <LeMonADE/Version.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/utility/ResultFormattingTools.h>

/***********************************************************************/
/**
 * @file
 *
 * @class AnalyzerWriteBfmFile
 *
 * @brief Analyzer writing the configurations into a given bfm-file.
 *
 * @details The output is appended to the file, if the file already exists.
 * If it does not exist, a new file is created and the header information is written
 * at the beginning
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @todo rename to WriteBFMFile or similar.
 **/
template <class IngredientsType>
class AnalyzerWriteBfmFile: public AbstractAnalyzer
{
public:
	//! Standard Constructor. Default: all data will be appended to the file
  AnalyzerWriteBfmFile(const std::string& filename,const IngredientsType& ing,int writeType=APPEND);

  //!Standard destructor closing the file stream.
  virtual ~AnalyzerWriteBfmFile();

  /**
   * @enum BFM_WRITE_TYPE
   *
   * @brief Specifiers for bfm-file-handling how data should be write-out.
   **/
  enum BFM_WRITE_TYPE{
	APPEND=0,	//!< The configuration (excl. header) is append to the file
	NEWFILE=1,	//!< The configuration (incl. header) is written to a new file
	OVERWRITE=2 //!< The configuration (incl. header) overwrites the existing file
  };
  
  /**
   * @enum BFM_WRITE_TYPE_EXPANDED
   *
   * @brief Specifiers the write type for the command handling of !add_bonds 
   * and !remove_bonds
   **/
  enum BFM_WRITE_TYPE_EXPANDED{
	  C_APPEND=3,	//!< The configuration (excl. header) is append to the file  
	  C_OVERWRITE=4,  //!< The configuration (incl. header) overwrites the existing file
	  C_NEWFILE=5,	//!< The configuration (incl. header) is written to a new file
	  C_APPNOFILE=6,	//!< The file doenst exist
  };
  
  

  //! Triggers writing of the Writes that need to be written for every step (i.e. excluding header)
  virtual bool execute();

  //! Initializes the writing. Opens the file and check for IO.
  virtual void initialize();

  //! Interface for adding new write Writes provided by Feature (e.g. !box_x)
  void registerWrite(std::string WriteString,SuperAbstractWrite* WriteObject);

  //! Interface for replacing old (standard) Writes with new write Writes. Example: FeatureJumps (!mcs)
  void replaceWrite(std::string WriteString,SuperAbstractWrite* NewWriteObject);

  //! Creates a new file (maybe overrides existing file)
  void startNewFile();

  //! Creates a new file with given \a name (maybe overrides existing file)
  void startNewFile(std::string fname);

  //! Creates a new file with given \a name and overrides existing file
  void startOverwriteNewFile(std::string fname);

  //! Closing the file stream.
  void closeFile(){file.close();}

  //! Returns a reference of Ingredients for general purpose
  const IngredientsType& getIngredients_() const {return ingredients;}

  //! Returns the filename used in this class
  std::string getFilename(){return _filename;}

  //! Returns the enum type for command handling
  int getCommandType(){return myCommandWriteType;}

  /**
   * @brief This function is called \a once in the end of the TaskManager.
   *
   * @details It´s a virtual function for inheritance.
   * Use this function for cleaning tasks (e.g. destroying arrays, file outut)
   *
   **/
  virtual void cleanup(){};
private:
  //! Write Writes that only need to be in the header
  void writeHeader();

  //! Checks if the file with given name already exists
  bool fileExists(std::string fname);

  //! The filename to be used.
  std::string _filename;

    //! Storage for data that are processed to file (mostly Ingredients).
  const IngredientsType& ingredients;

  //! The output file stream.
  std::ofstream file;

  //! ENUM-type BFM_WRITE_TYPE specify the write-out
  const int myWriteType;
  
  //! ENUM-type BFM_WRITE_TYPE_EXPANDED specify command handling
  int myCommandWriteType;

  //bfm Write strings with associated write objects
  //! Vector of pairs of Write-strings (e.g. !box_x) associated Write objects
  std::vector < std::pair<std::string, SuperAbstractWrite*> > WriteObjects;

  //init-flag to see if the initialize() routine was called. this
  //is necessary, because it opens the file and writes the header
  //! Flag for calling initialize() to open the file stream
  bool isInitialized;
};

/***********************************************************************/
//constructor
/***********************************************************************/
/**
 * @details It initialize all internal values and passes all information to the corresponding classes.
 * It register all known commands from the Feature.
 *
 *
 * @param filename Name of the file to write-out
 * @param ing Class holding all information of the system (mainly Ingredients )
 * @param writeType ENUM-type BFM_WRITE_TYPE to specify the write-out
 *
 * @todo Did I understand that correctly that the first !mcs is read not the last in the file?
 */
template <class IngredientsType>
AnalyzerWriteBfmFile<IngredientsType>::AnalyzerWriteBfmFile(const std::string& filename, const IngredientsType& ing, int writeType)
    :_filename(filename),ingredients(ing),myWriteType(writeType),isInitialized(false){}

/***********************************************************************/
//destructor
/***********************************************************************/
/**
 * @details it frees the memory from the write objects, which were registered
 * as pointers with registerWrite
 */
template<class IngredientsType>
AnalyzerWriteBfmFile<IngredientsType>::~AnalyzerWriteBfmFile()
{
	std::vector< std::pair <std::string,SuperAbstractWrite*> >::iterator it;

	for(it=WriteObjects.begin();it!=WriteObjects.end();++it)
	{
		delete it->second;
		it->second=0;
	}
	WriteObjects.clear();

	closeFile();

}

/***********************************************************************/
//void registerWrite
/***********************************************************************/
/**
 *
 * @param WriteString bfm-keyword/command/user-command (e.g. !box_x, #!fixed_monomers)
 * @param WriteObject Pointer to object providing the Write's functionality
 */
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::registerWrite(std::string WriteString,SuperAbstractWrite* WriteObject){


	std::vector< std::pair <std::string,SuperAbstractWrite*> >::iterator it;
	//first check if the command string (WriteString) is already used for
	//a different write
	for(it=WriteObjects.begin();it!=WriteObjects.end();++it)
	{
		if(it->first==WriteString)
		{

			std::stringstream errormessage;
			errormessage<<"WriteBfmFile::registerWrite "
					<<"Write string "<<WriteString<<" already used";
			throw std::runtime_error(errormessage.str());
		}
	}
	//if we are still here, we can register the write object
	WriteObjects.push_back(std::make_pair(WriteString,WriteObject));
	std::cout <<  WriteString << " registered for writing\n";
}


/***********************************************************************/
//void registerWrite
/***********************************************************************/
/**
 *
 * @param WriteString bfm-keyword/command/user-command (e.g. !box_x, #!fixed_monomers)
 * @param newWriteObject replace the write object with keyword key by the new one pointing to newWriteObject
 */
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::replaceWrite(std::string WriteString,SuperAbstractWrite* NewWriteObject)
{
	bool stringReplaced=false;

	std::vector< std::pair <std::string,SuperAbstractWrite*> >::iterator it;
	//first check if the command string (WriteString) is already used for
	//a different write
	for(it=WriteObjects.begin();it!=WriteObjects.end();++it)
	{
		if(it->first==WriteString)
		{
			//get the pointer to the old write object and free the memory
			SuperAbstractWrite* tmp= it->second;
			it->second=0;
			delete tmp;
			//insert the new object instead
			it->second=NewWriteObject;
			std::cout<<"WARNING: REPLACED WRITE "<<WriteString<<std::endl;
			stringReplaced=true;
		}
	}
	//if the write object was not replaced, throw an exception
	if(stringReplaced==false)
	{
		std::stringstream errormessage;
		errormessage<<"WriteBfmFile::replaceWrite "
			<<"Write string "<<WriteString<<" cannot be replaced, because it is not used yet!";
		throw std::runtime_error(errormessage.str());
	}

}

/***********************************************************************/
//boid startNewFile(string fname)
/***********************************************************************/
/**
 * @details Creates a new file with the existing filename inside this class.
 * It opens the file and writes the header. If the file already exists it adds
 * an "_" to the name to avoid overwriting of files.
 *
 * @throw <std::runtime_error> if general IO-error occurs
 **/
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::startNewFile()
{
    //this make a new file with the name _filename_1, _filename_2, etc.
    startNewFile(_filename);
}

/**
 * @details Creates a new file with the existing filename inside this class.
 * It opens the file and writes the header. If the file already exists it adds
 * an "_" to the name to avoid overwriting of files.
 *
 * @throw <std::runtime_error> if general IO-error occurs
 **/
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::startNewFile(std::string fname)
{
    //close the old file if necessary
    if(file.is_open()) file.close();
    //now determine the final name of the new file
    //if a file with this name exists, use fname_1 etc
    int i=1;
    std::string modifiedFileName(fname);

    while(fileExists(modifiedFileName) && i<10000)
    {
        std::stringstream manipulateName;
        manipulateName<<fname<<"_"<<i;
        ++i;
        modifiedFileName=manipulateName.str();
        std::cout<<modifiedFileName<<std::endl;
    }
    //set the _filename variable to the new value and open a new file.
    _filename=modifiedFileName;
    std::cout<<"opening new output file "<<_filename<<std::endl;
    file.open(_filename.c_str(),std::ios_base::in|std::ios_base::out|std::ios_base::app|std::ios_base::binary);
    if(file.fail()) throw std::runtime_error(std::string("WriteBfmFile: error opening output file ")+_filename);

    //write the bfm-header
    writeHeader();
}

/**
 * @details Creates a new file with the existing filename inside this class.
 * It opens the file and writes the header. If the file already exists it adds
 * an "_" to the name to avoid overwriting of files.
 *
 * @throw <std::runtime_error> if general IO-error occurs
 **/
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::startOverwriteNewFile(std::string fname)
{
    //close the old file if necessary
    if(file.is_open()) file.close();
    //now determine the final name of the new file
    //if a file with this name exists, use fname_1 etc

    std::cout<<"opening new output file "<<fname<<std::endl;
    file.open(fname.c_str(),std::ios_base::in|std::ios_base::out|std::ios_base::trunc|std::ios_base::binary);
    if(file.fail()) throw std::runtime_error(std::string("WriteBfmFile: error opening output file ")+fname);

    //write the bfm-header
    writeHeader();
}

/***********************************************************************/
//bool fileExists
/***********************************************************************/
/**
 * @param fname Filename for checking IO.
 * @return True if file already exists. False, otherwise.
 */
template <class IngredientsType>
bool AnalyzerWriteBfmFile<IngredientsType>::fileExists(std::string fname){
  //try to open the file for Writeing. if successful, the file exists
  std::ifstream f(fname.c_str());
  if(f.good()){
    f.close();
    return true;
  }
  else{
    f.close();
    return false;
  }
}

/***********************************************************************/
//bool execute
/***********************************************************************/
/**
 * @details Write-out of the next configuration of all Feature.
 *
 * @return True if everthing is alrigth. False if something goes wrong.
 */
template <class IngredientsType>
bool AnalyzerWriteBfmFile<IngredientsType>::execute(){

  if(myWriteType==OVERWRITE) startOverwriteNewFile(_filename);

  std::vector< std::pair<std::string,SuperAbstractWrite*> >::iterator it;
  SuperAbstractWrite* mcsCommand=0;
  //write all Writes that do not have the writeHeaderOnly flag set
  for(it=WriteObjects.begin(); it!=WriteObjects.end(); ++it)
  {
	if(it->first=="!mcs") mcsCommand=it->second;
	else if( (it->second)->writeHeaderOnly()==false)
	{
		(it->second)->writeStream(file);
	}
  }
  //if an !mcs was part of the write objects, write it now
  //the !mcs must always be written last in each step
  if(mcsCommand!=0) mcsCommand->writeStream(file);
  mcsCommand=0;
  return true;

}


/***********************************************************************/
//bool execute
/***********************************************************************/
/**
 * @details Handling the creation and openiong of file.
 * * If APPEND is specified and file exists all output is append to the existing file.
 * * If APPEND is specified and file don't exists, a new file is created and the header is written.
 * * If NEWFILE is specified and file don't exists don't exists, a new file is created and the header is written.
 * * If NEWFILE is specified and file exists, an exception is thrown.
 *
 * @throw <std::runtime_error> IO-error or Writing is not allowed or if unknown ENUM-type is used.
 */
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::initialize()
{
    if(myWriteType==APPEND && fileExists(_filename)) {myCommandWriteType=C_APPEND;}
    else if(myWriteType==APPEND && !fileExists(_filename)) {myCommandWriteType=C_APPNOFILE;}
    else if(myWriteType==NEWFILE) {myCommandWriteType=C_NEWFILE;}
    else if(myWriteType==OVERWRITE) {myCommandWriteType=C_OVERWRITE;}
    else/*flag incorrectly set*/
    {	std::cerr<<"AnalyzerWriteBfmFile:invalid flag " <<myCommandWriteType<<std::endl;
        throw std::runtime_error("WriteBfmFile: invalid flag set for writing. Valid options are APPEND or NEWFILE.\n");
    }
      //get the writing routines of all features
    ingredients.exportWrite(*this);
    
    //open the file. depending on whether or not the file exists,
    //and on the flag APPEND or NEWFILE a new file is created or not
    if(myWriteType==APPEND && fileExists(_filename))
    {
        std::cout<<"WriteBfmFile:appending to existing file "<<_filename<<std::endl;
        file.open(_filename.c_str(),std::ios_base::in|std::ios_base::out|std::ios_base::app|std::ios_base::binary);
        if(file.fail()) throw std::runtime_error(std::string("WriteBfmFile: error opening output file ")+_filename);
    }
    else if(myWriteType==APPEND && !fileExists(_filename))
    {
        std::cout<<"WriteBfmFile:file "<<_filename<<" does not exist. Creating new file instead"<<std::endl;
        startNewFile(_filename);
	
    }
    else if(myWriteType==NEWFILE)
    {
        if(!fileExists(_filename)) startNewFile(_filename);
        else
        {
            std::stringstream errormessage;
            errormessage<<"WriteBfmFile:trying to create new file "<<_filename
                       <<" on startup, but the file already exists."
                      <<"Choose a different filename or use APPPEND option in constructor of WriteBfmFile."<<std::endl;
            throw std::runtime_error(errormessage.str());
        }

    }
    else if(myWriteType==OVERWRITE)
    {
    	startOverwriteNewFile(_filename);

    }
    else/*flag incorrectly set*/
    {	
        throw std::runtime_error("WriteBfmFile: invalid flag set for writing. Valid options are APPEND or NEWFILE.\n");
    }
    isInitialized=true;
   
}

/***********************************************************************/
//void writeHeader
/***********************************************************************/
/**
 * @details It writes the Header of the file but not the first configuration.
 */
template <class IngredientsType>
void AnalyzerWriteBfmFile<IngredientsType>::writeHeader()
{
	file<<"########################################################\n";
	std::stringstream versionInfo;
	versionInfo<<"#!version="<<std::fixed<<std::setprecision(1)<<::LEMONADE_VERSION;
	file<<versionInfo.str()<<std::endl;
	file<<"#Bond Fluctuation Model - simulation data file"<<std::endl;
	time_t localTime=std::time(0);
	file<<"#"<<ctime(&localTime);
	file<<"#monomer numbering starts at 1\n\n";
	file<<"########################################################\n\n";

	file<<"# meta-data:"<<std::endl;
	std::stringstream metadata;
	ingredients.printMetaData(metadata);
	ResultFormattingTools::addComment(metadata);

	file<<metadata.str();
	file<<"########################################################\n\n";

	std::vector<std::pair<std::string,SuperAbstractWrite*> >::iterator it;
	//write all Writes, that have the writeHeaderOnly flag set
	for(it=WriteObjects.begin(); it!=WriteObjects.end(); ++it)
	{
		if( (it->second)->writeHeaderOnly())
		{
			(it->second)->writeStream(file);

		}
	}
}

#endif
