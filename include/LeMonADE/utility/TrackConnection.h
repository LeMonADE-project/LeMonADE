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

#ifndef LEMONADE_UPDATER_TRACKCONNECTION_H
#define LEMONADE_UPDATER_TRACKCONNECTION_H

#include <vector>
#include <sstream>
// #include <ostream>
#include <LeMonADE/utility/ResultFormattingTools.h>
/****************************************************************************/
template<class T=uint32_t >
class TrackLinks
{
public:
  
  TrackLinks():bufferSize(100),filename("BondCreationBreaking.dat"), isFirstFileDump(true), InformationSize(0) {};
//   TrackLinks(std::string filename_ = "BondCreationBreaking.dat" ){};
  
  /** 
   * @brief add reaction to be tracked
   * @details time in MCS and type is 0 for a break and a 1 for a new bond, monomer id's start at 0 
   */
  void addConnection(std::vector<T> vectorInformation) 
  { 
    for (uint32_t i = 0 ; i < InformationSize ; i++)
      Connection[i].push_back(vectorInformation[i]);
    if (Connection[0].size() == bufferSize) dumpReactions();
  }
  //! setter function for the size of the information to be stored 
  void setInformationSize(uint32_t InformationSize_)
  { 
    InformationSize=InformationSize_;
    resetConnection();
  } 
  
  //! getter function for the size of the information to be stored 
  uint32_t getInformationSize(){ return InformationSize;} 
  
  //! getter function for the connection size 
  const uint32_t getConnectionSize() const {return Connection.size();}

  //! dumps the data of the connection process into a file 
  void dumpReactions();  
  
  //! getter function for filename
  const std::string getFilename() const {return filename;}
  
  //! setter function for filename
  void setFilename(std::string filename_){filename=filename_;}
  
  //! getter function for the buffer size 
  const uint32_t getBufferSize() const {return bufferSize;}
  
  //! setter function for the buffer size
  void setBufferSize(uint32_t bufferSize_){bufferSize=bufferSize_;}
  
  //! add comment to the file header 
  void addComment(std::string comment_) 
  {
      if ( isFirstFileDump ) comment=comment_;
      else 
      {
	throw std::runtime_error("TrackLinks::addComment Try to add a comment, but the header is already written to the file!");
      }
  }
  
private:
  inline  bool resetConnection()
  {
    if (InformationSize > 0 )
    {
      Connection.resize(0);
      Connection.resize(InformationSize, std::vector<T>(T()) );
      return true;
    }else 
      throw std::runtime_error("TrackLinks::resetConnection size of the container is not set!");
      
  }
  //! container holding the information 
  std::vector< std::vector<T> > Connection;
  //!output filename 
  std::string filename;
  //! max length of internal buffer for each coordinate, before saving to disk
  uint32_t bufferSize;
  //! flag used in dumping time series output
  bool isFirstFileDump;
  //! comment for the header
  std::string comment;
  //! size(number of entries) of the information to be stored 
  uint32_t InformationSize;
};
/****************************************************************************/
template< class T >
void TrackLinks<T>::dumpReactions()
{
  	//fist make a single vector<vector<double> > for writing the results
	std::vector<std::vector<T> > resultsTimeseries=Connection;

	//if it is written for the first time, include comment in the output file
	if(isFirstFileDump){
	  
		std::ofstream file;
		file.open(filename.c_str());

		if(!file.is_open())
			throw std::runtime_error("TrackLinks::dumpReactions(): error opening output file"+filename+"\n");
		std::stringstream contents;
		ResultFormattingTools::writeTable(contents, resultsTimeseries, comment);

		file << contents.str();

		file.close();

		isFirstFileDump=false;
	}
	//otherwise just append the new data
	else{
		ResultFormattingTools::appendToResultFile(filename,
							  resultsTimeseries);
	}
	//set all time series vectors back to zero size
	resetConnection();
}
/****************************************************************************/

#endif 	/*LEMONADE_UPDATER_TRACKCONNECTION_H*/
