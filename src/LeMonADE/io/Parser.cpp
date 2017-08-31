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

#include <LeMonADE/io/Parser.h>

/*****************************************************************************/
/**
 * @file
 * @brief Implementation of Parser
 * */
/*****************************************************************************/

/*****************************************************************************/
//constructor and destructor
Parser::Parser(std::istream& inputStream)
  :stream(inputStream)
{}

Parser::~Parser(){}
/*****************************************************************************/

/*****************************************************************************/
//finds the next Read in the stream and returns the Readstring
std::string Parser::findRead()
{
	std::string line, Read;
	std::streampos linestart;

  while(!stream.eof() && !stream.fail()){
    //cout<<"findRead"<<endl;
    linestart=stream.tellg();
    getline(stream,line);
    if(!line.empty() && line.size()>1){
    bool ReadFound=(line.at(0)=='!' || (line.at(0)=='#' && line.at(1)=='!'));

    if(ReadFound){
      size_t length;
      //look for = sign
      length=line.find("=");
      if (length==std::string::npos){
	stream.seekg(linestart);
	stream>>Read;
	return Read;
      }
      else{
	stream.seekg(linestart);
	getline(stream,Read,'=');
	return Read;
      }

    }
    }
  }
  //if still here, the end of the file has been reached
  return "endoffile";
}
