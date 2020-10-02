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

#ifndef LEMONADE_UTILITY_RESULTFORMATTINGTOOLS_H
#define LEMONADE_UTILITY_RESULTFORMATTINGTOOLS_H

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <stdexcept>



/**
 * @file
 *
 * @namespace ResultFormattingTools
 *
 * @brief Helper function for formating and writing text output
 *
 * @details It provides IO-operation, writing all system information and the used Feature.
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @tparam ResultType Vector of results to write mainly( std::vector< std::vector <double> > ).
 *
 *
 * @todo we should reconsider this approach for usability
 * @todo move that to utility
 * @todo why we use a namespace?
 * @todo Example-Usage?
 * @todo Please comment...
 **/
namespace ResultFormattingTools {

/**
 * @brief Writes all given results into a given stream as table.
 *
 * @param stream Output stream
 * @param results List of results to format and write
 * @param comment Additional comments for the output
 */
template<class ResultType>
void writeTable(std::ostream& stream,
		ResultType results, std::string comment = "\n");

/**
 * @brief Adding a comment to the stream.
 * @param stream Output stream
 */
void addComment(std::stringstream& stream);

/**
 * @brief Writes all given results into a given stream as formatted output.
 *
 * @param filename Specify the name of the output-file.
 * @param ingredients A reference to the IngredientsType - mainly the system
 * @param results List of results to format and write
 * @param comment Additional comments for the output
 */
template<class IngredientsType, class ResultType>
void writeResultFile(std::string filename,const IngredientsType& ingredients,
		ResultType& results, std::string comment = "\n");

/**
 * @brief Appends all given results into a given stream as formatted output.
 *
 * @param filename Specify the name of the output-file.
 * @param results List of results to format and write
 */
template<class ResultType>
void appendToResultFile(std::string filename,ResultType& results);
}





template<class ResultType>
void ResultFormattingTools::writeTable(std::ostream& stream,
		ResultType results, std::string comment) {

	std::stringstream commentStream(comment);
	addComment(commentStream);
	stream << commentStream.str() << std::endl;

	//check if all column have the same size
	size_t columnSize = results[0].size();
	for (size_t i = 0; i < results.size(); ++i) {
		if (results[i].size() != columnSize)
			throw std::runtime_error("ResultFormattingTools::writeTable():Columns do not have the same size\n");
	}

	for (size_t row = 0; row < results[0].size(); ++row) {
		for (size_t column = 0; column < results.size(); ++column) {
			stream << results[column][row] << "\t";
		}
		stream << std::endl;
	}

}




template<class IngredientsType, class ResultType>
void ResultFormattingTools::writeResultFile(std::string filename,const IngredientsType& ingredients,
		ResultType& results, std::string comment) {

	std::ofstream file;
	file.open(filename.c_str());

	if(!file.is_open())
		throw std::runtime_error("ResultFormattingTools::writeResultFile(): error opening output file"+filename+"\n");
	std::stringstream contents;

	// write Header

	ingredients.printMetaData(contents);

	addComment(contents);

	writeTable(contents, results, comment);

	file << contents.str();

	file.close();
}

template<class ResultType>
void ResultFormattingTools::appendToResultFile(std::string filename,ResultType& results) {

	std::ofstream file;
	file.open(filename.c_str(),std::ios_base::app);

	if(!file.is_open())
		throw std::runtime_error("ResultFormattingTools::appendToResultFile(): error opening output file"+filename+"\n");

	//check if all column have the same size
	size_t columnSize = results[0].size();
	for (size_t i = 0; i < results.size(); ++i) {
		if (results[i].size() != columnSize){
			std::stringstream errormessage;
			errormessage<<"ResultFormattingTools::appendToResultFile():Columns do not have the same size\n";
			errormessage<<"first colums size "<<columnSize<<" column no "<<i<<" size "<<results[i].size()<<std::endl;
			throw std::runtime_error(errormessage.str());
		}
	}

	//write content
	for (size_t row = 0; row < results[0].size(); ++row) {
		for (size_t column = 0; column < results.size(); ++column) {
			file << results[column][row] << "\t";
		}
		file << std::endl;
	}

	file.close();
}

#endif /* LEMONADE_UTILITY_RESULTFORMATTINGTOOLS_H */
