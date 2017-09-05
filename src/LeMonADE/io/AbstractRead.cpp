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

#include <LeMonADE/io/AbstractRead.h>


/***********************************************************************
 * check if the string argument is a Read string (i.e ! or #! at beginning)
 ***********************************************************************/
/**
 * @brief Checks if given line contains a Read-string (i.e. !... or #!...)
 *
 * @details This function tests the given line if the first characters are an \a !
 * \a #! representing a standard command or a user defined command.
 *
 * @param line A string to detect if command are present.
 * @return True if line is an command. False - everything else.
 */
bool AbstractRead::detectRead(std::string& line) const {
    if (line.at(0)=='!')
        return true;
    else if (line.size()>1 && line.at(0)=='#' && line.at(1)=='!')
        return true;
    else
        return false;
}

/***********************************************************************
 * checks if the next character in the stream is separator. ignores whitespace.
 ***********************************************************************/
/**
 * @brief Checks if the next character in the stream is separator ignoring whitespace.
 *
 * @param stream Input stream to operate on.
 * @param separator Look for this character.
 *
 * @return True if the next character is the seperator. False otherwise.
 */
bool AbstractRead::findSeparator(std::istream& stream, char separator) {
    while (stream.peek()==' ' && stream.good()) {
        stream.ignore(1);
    }
    if (stream.peek()==separator) {
        stream.ignore(1);
        return true;
    }
    else
        return false;
}

/******************************************************************************
 * tokenizes the string str into three strings. the separators are given by delim1,2
 * ***************************************************************************/
/**
 * @brief This function tokenize a string into tokens by 2 given delimiters.
 *
 * @param str String to be tokenized
 * @param delim1 first delimiter
 * @param delim2 second delimiter
 * @return vector<string> of tokens
 */
/***********************************************************************/
std::vector<std::string> AbstractRead::tokenize2Parameter(const std::string& str,
  		char delim1, char delim2) {
  	std::vector<std::string> tokens;
  	std::stringstream mySstream(str);
  	std::string temp;

  	while (getline(mySstream, temp, delim1)) {

  		if(temp.size() != 0)
  		{
  			std::stringstream mySstream2(temp);
  			std::string temp2;
  			while (getline(mySstream2, temp2, delim2)) {
  				//printf("%s \n", temp2.c_str());
  				if(temp2.size() != 0)
  					tokens.push_back(temp2);
  			}
  		}

  	}

  	return tokens;

  }
