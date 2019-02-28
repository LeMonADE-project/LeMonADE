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
#include <ostream>
/****************************************************************************/
struct reaction
{
  reaction(uint32_t time_, bool type_, uint32_t MonID1_, uint32_t MonID2_):time(time_),type(type_),MonID1(MonID1_), MonID2(MonID2_){};
  uint32_t time;
  bool type;
  uint32_t MonID1;
  uint32_t MonID2;
  /****************************************************************************/
  /**
  * @brief \b Stream \b Out \b operator of the reaction
  *
  * @details Streams out the elements of the class 
  *
  * @param stream output-stream
  * @param label object of class reaction
  * @return output-stream
  **/
  friend std::ostream& operator<< (std::ostream& stream, const reaction & Reactivity)
  {
	  stream 	  
	  << Reactivity::type   << " " 
	  << Reactivity::time   << " "
	  << Reactivity::MonID1 << " "
	  << Reactivity::MonID2 << " ";
	  return stream;
  }
};

/****************************************************************************/
class tracker
{
public:
  
  tracker():filename("BondCreationBreaking.dat"){};
  tracker(std::string filename_ = "BondCreationBreaking.dat" ){};
  
  /** 
   * @brief add reaction to be tracked
   * @details time in MCS and type is 0 for a break and a 1 for a new bond, monomer id's start at 0 
   */
  addReaction(uint32_t time_, bool type_, uint32_t MonID1_, uint32_t MonID2_)
  {
    ConnectionsBreaks.push_back(reaction( time_, type_, MonID1_, MonID2_));
  }
  void setFilename(std::string filename_){filename=filename_;};
  std::string getFilename() {return filename;};
  
private:
  //! dumps the data of the connection process into a file 
  void dumpReactions();  
  //! container holding the information 
  std::vector< reaction > ConnectionsBreaks;
  //output filename 
  std::string filename;
};
/****************************************************************************/
void tracker::dumpReactions()
{
  std::ofstream reactbreakFile;
  reactbreakFile.open(filename.c_str(),std::ios_base::app);
  for(size_t n=0;n<ConnectionsBreaks.size();n++)
    reactbreakFile<<ConnectionsBreaks[n]<<std::endl;
  reactbreakFile.close();
  ConnectionsBreaks.clear();
}
/****************************************************************************/

#endif 	/*LEMONADE_UPDATER_TRACKCONNECTION_H*/
